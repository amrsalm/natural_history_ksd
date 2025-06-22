# Life-Course Comorbidity Analysis in Kidney Stone Disease (KSD)

This repository contains code and visualizations for a phenome-wide, age-stratified analysis of comorbidity burden in individuals with kidney stone disease (KSD), using electronic health record (EHR) data. The methods focus on both **odds ratio-based enrichment analysis** and **absolute risk difference (ARD) visualizations** to uncover systemic risk patterns across the lifespan.

---

## Repository Structure

analysis/
│
├── data/
│ └── results.csv # Cleaned input dataset containing diagnosis counts
│
├── odds_ratio_analysis.ipynb # Notebook for computing odds ratios per ICD-10 chapter
│ # across age intervals, with FDR correction and heatmap
│
├── ard_bubble_plot.py # Script for generating ARD-based bubble plots of top
│ # phenotypes (e.g., neoplasms), stratified by age
│
├── figures/
│ ├── figure1_odds_heatmap.pdf # Output: heatmap of chapter-level odds ratios
│ └── figure2_neoplasm_bubbles.pdf # Output: ARD bubble plot for neoplasms (ages 70–80)
│
└── utils/
└── preprocessing.py # Data wrangling utilities: filtering, code propagation,
# and binning across age intervals

---

## Key Analyses

### 1. Odds Ratio Analysis

- **Purpose**: Quantify overrepresentation of ICD-10 chapters in KSD vs matched controls
- **Method**:
  - Performed two-tailed Fisher's exact tests across 10-year age intervals
  - Applied Benjamini–Hochberg FDR correction
- **Output**: Heatmap of statistically significant associations (Figure 1)

---

### 2. Absolute Risk Difference (ARD) Bubble Plot

- **Purpose**: Visualize the most common phenotypes within a disease domain (e.g., neoplasms)
- **Method**:
  - Compute ARD:

    \[
    \text{ARD} = \frac{n_{\text{stone}}}{N_{\text{stone}}} - \frac{n_{\text{no stone}}}{N_{\text{no stone}}}
    \]

  - Bubble size represents the number of stone formers with that phenotype
  - Marker shape distinguishes assigned vs propagated diagnoses
- **Output**: Bubble plot (Figure 2) for top neoplastic phenotypes in ages 70–80

---

## Requirements

- Python ≥ 3.8
- pandas
- seaborn
- matplotlib
- scipy
- statsmodels

Install dependencies using:


