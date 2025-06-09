import pandas as pd

import pandas as pd
import ast
import pandas as pd
import numpy as np
import scipy.stats as stats
import ast


merged_CASE_df = pd.read_csv("casetotest_0912.csv")
merged_Control_Df = pd.read_csv("controltotest_0912.csv")
ICD_prop = pd.read_csv("icd10_propagation (1) (1).csv")


# Convert the 'ancestors' column into a list (array)
ICD_prop['ancestors'] = ICD_prop['ancestors'].str.split(',')

merged_CASE_df = pd.merge(merged_CASE_df, ICD_prop[['code', 'ancestors']], left_on='diagnosis_code', right_on='code', how='left')
merged_Control_Df =pd.merge(merged_Control_Df, ICD_prop[['code', 'ancestors']], left_on='diagnosis_code', right_on='code', how='left')
# Now, we add the 'ancestors' to the 'parents' column, handling the case where 'parents' is NaN (root nodes)
merged_Case_Df_cleaned = merged_CASE_df.iloc[:, 1:]
merged_Control_Df_cleaned = merged_Control_Df.iloc[:, 1:]
# Remove rows where 'parents' column is NaN
merged_Case_Df_cleaned = merged_CASE_df.dropna(subset=['ancestors'])
merged_Control_Df_cleaned = merged_Control_Df.dropna(subset=['ancestors'])

icd10_chapter_mapping = [
    ('A00', 'B99', 1),  # Certain infectious and parasitic diseases
    ('C00', 'D49', 2),  # Neoplasms
    ('D50', 'D89', 3),  # Diseases of the blood and blood-forming organs
    ('E00', 'E89', 4),  # Endocrine, nutritional and metabolic diseases
    ('F01', 'F99', 5),  # Mental, Behavioral and Neurodevelopmental disorders
    ('G00', 'G99', 6),  # Diseases of the nervous system
    ('H00', 'H59', 7),  # Diseases of the eye and adnexa
    ('H60', 'H95', 8),  # Diseases of the ear and mastoid process
    ('I00', 'I99', 9),  # Diseases of the circulatory system
    ('J00', 'J99', 10), # Diseases of the respiratory system
    ('K00', 'K95', 11), # Diseases of the digestive system
    ('L00', 'L99', 12), # Diseases of the skin and subcutaneous tissue
    ('M00', 'M99', 13), # Diseases of the musculoskeletal system
    ('N00', 'N99', 14), # Diseases of the genitourinary system
    ('O00', 'O9A', 15), # Pregnancy, childbirth, and the puerperium
    ('P00', 'P96', 16), # Certain conditions originating in the perinatal period
    ('Q00', 'Q99', 17), # Congenital malformations
    ('R00', 'R99', 18), # Symptoms and abnormal findings
    ('S00', 'T88', 19), # Injury, poisoning
    ('U00', 'U85', 20), # Codes for special purposes
    ('V00', 'Y99', 21), # External causes of morbidity
    ('Z00', 'Z99', 22)  # Factors influencing health status
]
def get_chapter(row):
    # Choose 'Original_Phenotype' if available, otherwise 'Phenotype'
    code = row['diagnosis_code']

    if pd.notna(code):
        # Extract the first letter and first two digits (e.g., 'A00')
        code_prefix = code[:3].upper()

        # Match the code to the correct ICD-10 chapter range
        for start, end, chapter in icd10_chapter_mapping:
            if start <= code_prefix <= end:
                return chapter

    return np.nan  # Return NaN if no match is found or code is NaN
merged_Control_Df_cleaned['chapter'] = merged_Control_Df_cleaned.apply(get_chapter, axis=1)
merged_Case_Df_cleaned['chapter'] = merged_Case_Df_cleaned.apply(get_chapter, axis=1)
import pandas as pd
import numpy as np
import ast
from scipy import stats
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
from tqdm import tqdm

# Define Fisher's exact test
def calculate_fisher_exact(stone_count, nostone_count, total_stone_count, total_nostone_count):
    stone_no_phenotype = total_stone_count - stone_count
    nostone_no_phenotype = total_nostone_count - nostone_count
    contingency_table = [[stone_count, nostone_count], [stone_no_phenotype, nostone_no_phenotype]]
    _, p_value = stats.fisher_exact(contingency_table)
    return p_value



# Define age bins and labels
bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, np.inf]
bin_labels = ['0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '>90']

# Get EMR coverage
stone_coverage = merged_Case_Df_cleaned.groupby('primary_mrn')['diagnosis_age'].agg(['min', 'max']).reset_index()
stone_coverage.columns = ['primary_mrn', 'first_encounter', 'last_encounter']
nostone_coverage = merged_Control_Df_cleaned.groupby('primary_mrn')['diagnosis_age'].agg(['min', 'max']).reset_index()
nostone_coverage.columns = ['primary_mrn', 'first_encounter', 'last_encounter']

# Create dict to store results
results = {'Interval': [], 'Phenotype': [], 'Original_Phenotype': [], 'P-Value': [],
           'Stone_Phentoype': [], 'NoStone_Phenotype': [], 'Stone_NoPhenotype': [],
           'NoStone_NoPhenotype': [], 'Total_Stone': [], 'Total_NoStone': []}

# Function for parallel row expansion
def expand_phenotypes(row):
    expanded = []
    mrn = row['primary_mrn']
    diagnosis_code = row['diagnosis_code']
    diagnosis_age = row['diagnosis_age']
    expanded.append((mrn, diagnosis_code, diagnosis_code, diagnosis_age))
    try:
        ancestors = ast.literal_eval(str(row.get('ancestors', '')))
        if isinstance(ancestors, list):
            for ancestor_set in ancestors:
                for ancestor in ancestor_set.split(';'):
                    expanded.append((mrn, ancestor, diagnosis_code, diagnosis_age))
    except (ValueError, SyntaxError):
        pass
    return expanded

# Expand with parallelism
def expand_all(df, desc="Expanding"):
    results = Parallel(n_jobs=8)(delayed(expand_phenotypes)(row) for _, row in tqdm(df.iterrows(), total=len(df), desc=desc))
    flat = [item for sublist in results for item in sublist]
    return pd.DataFrame(flat, columns=['primary_mrn', 'phenotype', 'original_phenotype', 'diagnosis_age'])

# Process by interval
for i in range(len(bins) - 1):
    interval_min, interval_max = bins[i], bins[i+1]
    interval_label = bin_labels[i]

    # Get in-system patients
    in_system_stone = stone_coverage[
        (stone_coverage['first_encounter'] <= interval_max) & (stone_coverage['last_encounter'] >= interval_min)
    ]
    in_system_nostone = nostone_coverage[
        (nostone_coverage['first_encounter'] <= interval_max) & (nostone_coverage['last_encounter'] >= interval_min)
    ]

    # Filter main dataframes
    df_stone_interval = merged_Case_Df_cleaned[
        (merged_Case_Df_cleaned['diagnosis_age'] >= interval_min) &
        (merged_Case_Df_cleaned['diagnosis_age'] < interval_max) &
        (merged_Case_Df_cleaned['primary_mrn'].isin(in_system_stone['primary_mrn']))
    ]
    df_nostone_interval = merged_Control_Df_cleaned[
        (merged_Control_Df_cleaned['diagnosis_age'] >= interval_min) &
        (merged_Control_Df_cleaned['diagnosis_age'] < interval_max) &
        (merged_Control_Df_cleaned['primary_mrn'].isin(in_system_nostone['primary_mrn']))
    ]

    # Expand phenotypes
    stone_pheno = expand_all(df_stone_interval, desc=f"Stone {interval_label}")
    nostone_pheno = expand_all(df_nostone_interval, desc=f"NoStone {interval_label}")

    # Unique patients
    total_stone = in_system_stone['primary_mrn'].nunique()
    total_nostone = in_system_nostone['primary_mrn'].nunique()

    # Phenotype universe
    all_phens = set(stone_pheno['phenotype']).union(nostone_pheno['phenotype'])

    for pheno in all_phens:
        stone_count = stone_pheno[stone_pheno['phenotype'] == pheno]['primary_mrn'].nunique()
        nostone_count = nostone_pheno[nostone_pheno['phenotype'] == pheno]['primary_mrn'].nunique()
        pval = calculate_fisher_exact(stone_count, nostone_count, total_stone, total_nostone)

        original = stone_pheno[stone_pheno['phenotype'] == pheno]['original_phenotype'].iloc[0] \
            if not stone_pheno[stone_pheno['phenotype'] == pheno].empty else None

        results['Interval'].append(interval_label)
        results['Phenotype'].append(pheno)
        results['Original_Phenotype'].append(original)
        results['P-Value'].append(pval)
        results['Stone_Phentoype'].append(stone_count)
        results['NoStone_Phenotype'].append(nostone_count)
        results['Stone_NoPhenotype'].append(total_stone - stone_count)
        results['NoStone_NoPhenotype'].append(total_nostone - nostone_count)
        results['Total_Stone'].append(total_stone)
        results['Total_NoStone'].append(total_nostone)

# Create DataFrame and adjust p-values
df_results = pd.DataFrame(results)
df_results['P-Value_Adjusted'] = multipletests(df_results['P-Value'], method='fdr_bh')[1]
df_results['Interval'] = df_results['Interval'].replace(">90", "90-120")

# Save
df_results.to_csv("results_in_system_parallel.csv", index=False)
