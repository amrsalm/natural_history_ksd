import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Assume your data is called most_significant_associations

# 2. Target ICD-10 codes
target_codes = [
    "E70-E88",   # Metabolic disorders
    'I10-I1A',   # Hypertensive diseases
    'N17-N19',   # Acute kidney failure and CKD
    'E08-E13',   # Diabetes mellitus
    'E65-E68',   # Obesity and hyperalimentation
'19',
    'Z72.0'

]

# 3. Make sure Phenotype column is string and cleaned (strip spaces!)
most_significant_associations['Phenotype'] = most_significant_associations['Phenotype'].astype(str).str.strip()

# 4. Subset rows where Phenotype matches exactly any of the target codes
df_target = most_significant_associations[
    most_significant_associations['Phenotype'].isin(target_codes)
].copy()


# --- Extract starting age from Interval to sort numerically ---
# This assumes intervals are in the form like "10-20", "20-30", ">90", etc.
def extract_start_age(interval):
    if isinstance(interval, str):
        if interval.startswith('>'):
            return int(interval[1:]) + 1
        try:
            return int(interval.split('-')[0])
        except:
            return np.nan
    return np.nan

# Apply to create sort key
most_significant_associations['Interval_Start'] = most_significant_associations['Interval'].apply(extract_start_age)

# Clean and prepare again with sorting
df_target = most_significant_associations[
    most_significant_associations['Phenotype'].isin(target_codes)
].copy()

df_target['Risk_Difference'] = (
    (df_target['Stone_Phentoype'] / (df_target['Total_Stone'] )) -
    (df_target['NoStone_Phenotype'] / (df_target['Total_NoStone'] ))
)

df_target = df_target.replace([np.inf, -np.inf], np.nan)
df_target = df_target.dropna(subset=['Risk_Difference'])

# Sort first by Interval_Start, then by Phenotype
df_target = df_target.sort_values(['Interval_Start', 'Phenotype'])

# Recreate mapping
phenotype_label_map = df_target.drop_duplicates('Phenotype').set_index('Phenotype')['code_description'].to_dict()

# Plot
plt.figure(figsize=(14, 8))

# Get sorted unique intervals for consistent x-axis ordering
ordered_intervals = df_target[['Interval', 'Interval_Start']].drop_duplicates().sort_values('Interval_Start')['Interval'].tolist()

# Plot setup

for code in target_codes:
    subset = df_target[df_target['Phenotype'] == code]
    if not subset.empty:
        label = phenotype_label_map.get(code, code)

        # Reindex to ensure consistent x-axis
        subset = subset.set_index('Interval').reindex(ordered_intervals).reset_index()

        # Plot line first (connect all points)
        plt.plot(subset['Interval'], subset['Risk_Difference'], label=label, linestyle='-', linewidth=2)

        # Plot each point individually based on significance
        for _, row in subset.iterrows():
            if pd.isna(row['Risk_Difference']):
                continue

            marker_style = 'o' if row['P-Value_Adjusted'] < 0.05 else 'x'
            marker_color = 'black'
            plt.plot(row['Interval'], row['Risk_Difference'], marker=marker_style, color=marker_color)


# Decorate
plt.xlabel('Age Interval', fontsize=14)
plt.ylabel('Risk Difference (Stone vs. Control)', fontsize=14)
plt.title('Absolute Risk Difference of Disease in Stone Formers by Age', fontsize=16)
plt.legend(title='Disease', fontsize=11, title_fontsize=12)
plt.xticks(rotation=45)
plt.grid(True)
plt.tight_layout()
# Custom legend handles
import matplotlib.lines as mlines

significant_marker = mlines.Line2D([], [], color='black', marker='o', linestyle='None', label='Significant (p < 0.05)')
nonsignificant_marker = mlines.Line2D([], [], color='gray', marker='x', linestyle='None', label='Not Significant')

plt.legend(handles=[*plt.gca().get_legend_handles_labels()[0], significant_marker, nonsignificant_marker],
           title='Disease / Significance', fontsize=11, title_fontsize=12)

# Save
plt.savefig('risk_difference_by_age_sorted_withchatper19_marker.png', dpi=300)
