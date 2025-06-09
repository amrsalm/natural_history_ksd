import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np
import plotly.graph_objs as go

# Assuming `most_significant_associations` and `df_significant` are your DataFrames

# Convert to string first if not already
most_significant_associations['Phenotype'] = most_significant_associations['Phenotype'].astype(str)



# Sort by 'chapter' to ensure chapters are sorted numerically on the y-axis
most_significant_associations = most_significant_associations.sort_values('chapter')

chapter_data = most_significant_associations
chapter_data= chapter_data.dropna(subset=['chapter'])
# Define the ICD-10 chapter names
icd_chapter_names = {
    1: "Infectious Diseases",
    2: "Neoplasms",
    3: "Blood Disorders",
    4: "Endocrine Diseases",
    5: "Mental Disorders",
    6: "Nervous System",
    7: "Eye Diseases",
    8: "Ear Diseases",
    9: "Circulatory Diseases",
    10: "Respiratory Diseases",
    11: "Digestive Diseases",
    12: "Skin Diseases",
    13: "Musculoskeletal Diseases",
    14: "Genitourinary Diseases",
    15: "Pregnancy Complications",
    16: "Perinatal Conditions",
    17: "Congenital Abnormalities",
    18: "Symptoms and Signs",
    19: "Injury and Poisoning",
    20: "External Causes",
    21: "Health Status",
    22: "Special Purposes"
}

# Add a new column 'chapter_name' with the corresponding names
chapter_data['chapter_name'] = chapter_data['chapter'].map(icd_chapter_names)


# Step 3: Compute the odds ratio
chapter_data['Odds_Ratio_Computed'] = (chapter_data['Stone_Phentoype'] / chapter_data['NoStone_Phenotype']) / \
                                      (chapter_data['Stone_NoPhenotype'] / chapter_data['NoStone_NoPhenotype'])

# Step 4: Handle infinite values by setting them to a very large number
chapter_data['Odds_Ratio_Computed'].replace([np.inf, -np.inf], np.nan, inplace=True)

chapter_data_phenotype = chapter_data.copy()
# Filter rows where 'Phenotype' is exactly between 1 and 22
chapter_data = chapter_data[
    chapter_data['Phenotype'].isin([str(i) for i in range(1, 23)])
]

# Step 2: Convert the 'Phenotype' and 'chapter' columns to integers for matching
chapter_data['Phenotype'] = chapter_data['Phenotype'].astype(int)
chapter_data['chapter'] = chapter_data['chapter'].astype(int)
chapter_data_phenotype['chapter'] = chapter_data_phenotype['chapter'].astype(int)
# Step 5: Calculate percentiles (25th, 50th, 75th) for finite values
q25 = np.percentile(chapter_data['Odds_Ratio_Computed'].dropna(), 25)
median = np.percentile(chapter_data['Odds_Ratio_Computed'].dropna(), 50)
q75 = np.percentile(chapter_data['Odds_Ratio_Computed'].dropna(), 75)
finite_min = chapter_data['Odds_Ratio_Computed'].dropna().min()
finite_max = chapter_data['Odds_Ratio_Computed'].dropna().max()

# Step 6: Replace 'inf' with a value slightly larger than the maximum finite value for coloring purposes
inf_value = finite_max * 1.5
chapter_data['Odds_Ratio_Computed'].replace(np.nan, inf_value, inplace=True)
# Log transform odds ratios
log_odds = np.log1p(heatmap_data)

# Set up color scale
vmin = np.nanmin(log_odds.values[np.isfinite(log_odds.values)])
vmax = np.nanmax(log_odds.values[np.isfinite(log_odds.values)])
center = np.log1p(1)  # corresponds to OR=1

norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
cmap = mcolors.LinearSegmentedColormap.from_list("custom_bwr", ["dodgerblue", "white", "firebrick"])

# Create figure
fig, ax = plt.subplots(figsize=(18, 10))

# Draw colored boxes manually
for i, chapter in enumerate(heatmap_data.index):
    for j, interval in enumerate(heatmap_data.columns):
        val = heatmap_data.loc[chapter, interval]

        if pd.isna(val):
            color = "lightgray"
            log_val = np.nan
        elif val == np.inf:
            # Use highest color in colorW scale for inf
            log_val = vmax
            color = cmap(norm(vmax))
        else:
            log_val = np.log1p(val)
            color = cmap(norm(log_val))


        is_sig = ((df_significant['chapter'] == chapter) & (df_significant['Interval'] == interval)).any()

        # Draw rectangle with or without hatch based on significance
        hatch = None if is_sig else '////'
        edgecolor = 'white' if not is_sig else 'white'
        rect = mpatches.Rectangle((j, i), 1, 1, facecolor=color, edgecolor=edgecolor,
                                  hatch=hatch, linewidth=0.5)
        ax.add_patch(rect)

        # Determine annotation text
        if pd.isna(val):
            text = "N/A"
            font_color = "black"
        elif val == inf_value:
            text = "∞" if is_sig else "∞"
            font_color = "black"
        else:
            stars = "**" if is_sig else ""
            text = f"{val:.2f}{stars if not is_sig else ''}"
            font_color = "black" if val >= 1 else "blue"

        ax.text(j + 0.5, i + 0.5, text, ha='center', va='center', fontsize=10, color=font_color)

# Set ticks and labels
ax.set_xlim(0, len(heatmap_data.columns))
ax.set_ylim(0, len(heatmap_data.index))

ax.set_xticks(np.arange(len(heatmap_data.columns)) + 0.5)
ax.set_xticklabels(heatmap_data.columns, rotation=45, ha='right')

ax.set_yticks(np.arange(len(heatmap_data.index)) + 0.5)
ax.set_yticklabels([icd_chapter_names[i] for i in heatmap_data.index])

ax.invert_yaxis()
ax.set_xlabel("Age Interval")
ax.set_ylabel("ICD-10 Chapter")
ax.set_title("Odds Ratio Heatmap", fontsize=16)

# Create colorbar manually
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # dummy array for colorbar
cbar = plt.colorbar(sm, ax=ax, orientation='vertical', label='Odds Ratio')

plt.tight_layout()
plt.savefig("odds_ratio_heatmap_solid_for_significance.png", dpi=300)
plt.show()
