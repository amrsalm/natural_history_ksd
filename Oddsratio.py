import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# Step 8: Define the smoother color scale
colorscale = [
    [0, 'rgb(30, 144, 255)'],  # Sky Blue: rgb(135, 206, 235)
    [1, 'white'],              # White at the center
    [4, 'rgb(255, 111, 97)']     # Soft Red: rgb(255, 111, 97)
]
# Log transform odds ratios
log_odds = np.log(heatmap_data)

# Set up color scale
vmin = np.nanmin(log_odds.values[np.isfinite(log_odds.values)])
vmax = np.nanmax(log_odds.values[np.isfinite(log_odds.values)])
center = 0  # corresponds to OR=1

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
            log_val = np.log(val)
            color = cmap(norm(log_val))


        is_sig = ((df_significant_numeric['chapter'] == chapter) & (df_significant_numeric['Interval'] == interval)).any()

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
        elif val == inf_value :
            text = "∞"
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
log_ticks = np.array([np.log(0.25), np.log(0.5), np.log(1), np.log(2), np.log(4)])
tick_labels = ["0.25", "0.5", "1", "2", "4"]

cbar.set_ticks(log_ticks)
cbar.set_ticklabels(tick_labels)

plt.tight_layout()
plt.savefig("odds_ratio_heatmap_solid_for_significance.pdf", dpi=1000, format='pdf')
plt.show()
