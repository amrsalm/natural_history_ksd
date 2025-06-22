import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap

# Filter the data for neoplasms (chapter 2) and a specific interval
df_plot = results[
    (results['chapter'] == 2) &
    (results['Interval'] == '70-80')
].copy()

df_plot['is_propagated'] = df_plot['Phenotype'] != df_plot['Original_Phenotype']


# Calculate ARD and prevalence
df_plot['risk_stone'] = df_plot['Stone_Phentoype'] / df_plot['Total_Stone']
df_plot['risk_nostone'] = df_plot['NoStone_Phenotype'] / df_plot['Total_NoStone']
df_plot['ARD'] = df_plot['risk_stone'] - df_plot['risk_nostone']
df_plot['Prevalence'] = df_plot['Stone_Phentoype']

# Sort by prevalence and select top 10
df_plot = df_plot[df_plot['Prevalence'] > 0]
df_plot = df_plot.sort_values(by='Prevalence', ascending=False).head(10)

# Wrap long phenotype descriptions for y-axis labels
df_plot['code_description_wrapped'] = df_plot['code_description'].apply(
    lambda x: '\n'.join(textwrap.wrap(x.strip(), width=40))
)

# Plot settings
plt.figure(figsize=(14, 10))  # more vertical space
sns.set(style="whitegrid")

# Bubble plot
sns.scatterplot(
    data=df_plot,
    x='ARD',
    y='code_description_wrapped',
    size='Prevalence',
    hue='ARD',
    style='is_propagated',  # This automatically picks two marker styles
    palette=sns.diverging_palette(220, 20, as_cmap=True),
    sizes=(100, 1000),
    legend=False
)

# Add counts as text
for _, row in df_plot.iterrows():
    weight = 'bold' if row['is_propagated'] else 'normal'
    plt.text(row['ARD'] + 0.01, row['code_description_wrapped'], f"{row['Stone_Phentoype']}",
             fontsize=8, va='center', ha='left', color='black', fontweight=weight)


# Decorations
# Decorations
plt.axvline(0, color='gray', linestyle='--')
plt.xlabel('Absolute Risk Difference (ARD)', fontsize=12)
plt.ylabel('Neoplastic Phenotype', fontsize=12)
plt.title('Figure 2: Top 10 Neoplastic Phenotypes in Stone Formers (Age 70–80)', fontsize=16)

# Add vertical padding
plt.ylim(plt.ylim()[0] + 0.5, plt.ylim()[1] - 0.5)


plt.tight_layout(pad=3)
plt.savefig("figure2_neoplasm_top10_bubble_plot.pdf", dpi=1000, format='pdf')
plt.show()

