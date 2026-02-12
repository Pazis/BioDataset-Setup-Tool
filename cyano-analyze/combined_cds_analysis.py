import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

# --- CONFIGURATION ---
# List your CSVs here. Format: ("Path/To/File.csv", "FamilyName")
INPUT_FILES = [
    ("/home/angelos/Desktop/Thesis/test/full-mode/microcystaceae_summary.csv", "Microcystaceae"),
    ("/home/angelos/Desktop/Thesis/test/full-mode/nostocaceae_summary.csv", "Nostocaceae"),
    ("/home/angelos/Desktop/Thesis/test/full-mode/prochlorococaceae_summary.csv", "Prochlorococaceae"),
    ("/home/angelos/Desktop/Thesis/test/full-mode/synechococaceae_summary.csv", "Synechococcaceae")
]
OUTPUT_DIR = "combined_plots"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# --- 1. Load and Combine Data ---
all_data = []

def clean_pct(val):
    if isinstance(val, str):
        return float(val.replace('%', ''))
    return float(val)

for file_path, family in INPUT_FILES:
    if os.path.exists(file_path):
        temp_df = pd.read_csv(file_path)
        temp_df['Family'] = family  # Add a column to identify the group
        all_data.append(temp_df)
    else:
        print(f"Warning: Skipping {file_path}, file not found.")

if not all_data:
    print("Error: No data loaded. Check your file paths.")
    sys.exit(1)

df = pd.concat(all_data, ignore_index=True)
df['Initial_Undesc_PCT'] = df['Initial_Undesc_PCT'].apply(clean_pct)
df['Enhanced_Undesc_PCT'] = df['Enhanced_Undesc_PCT'].apply(clean_pct)
df['Rel_Reduction_PCT'] = (
    (df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']) / 
    df['Initial_Undesc_PCT'].replace(0, 1)
) * 100

# --- ΝΕΟΙ ΥΠΟΛΟΓΙΣΜΟΙ ---

# 1. Υπολογισμός Συνολικών Γονιδίων (Total Genes)
# Βασιζόμαστε στα Initial Hypotheticals και το Initial Percentage
# Χρησιμοποιούμε replace(0, 1) για αποφυγή διαίρεσης με το μηδέν
df['Total_Genes_Est'] = (df['Initial_Hypotheticals'] / (df['Initial_Undesc_PCT'] / 100)).replace(0, 1)

# 2. Υπολογισμός Ποσοστού Ψευδογονιδίων (%)
df['Initial_Pseudo_PCT'] = (df['Initial_Pseudogene_candidates'] / df['Total_Genes_Est']) * 100
df['Enhanced_Pseudo_PCT'] = (df['Enhanced_Pseudogene_candidates'] / df['Total_Genes_Est']) * 100

# Σχετική μείωση Undescribed (Impact)
df['Rel_Reduction_PCT'] = (
    (df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']) / 
    df['Initial_Undesc_PCT'].replace(0, 1)
) * 100
# ---------------------------------------------------------
# A. STATISTICAL ANALYSIS (Per Family)
# ---------------------------------------------------------
print("\n--- 1. Statistical Analysis & Impact Tracking ---")

for family in df['Family'].unique():
    fam_df = df[df['Family'] == family]
    
    avg_initial = fam_df['Initial_Undesc_PCT'].mean()
    avg_enhanced = fam_df['Enhanced_Undesc_PCT'].mean()
    
    # This is the metric you requested:
    total_rel_improvement = fam_df['Rel_Reduction_PCT'].mean()
    
    print(f"FAMILY: {family}")
    print(f"  - Initial Hypothetical:  {avg_initial:.2f}%")
    print(f"  - Enhanced Hypothetical: {avg_enhanced:.2f}%")
    print(f"  - IMPACT: Your pipeline resolved {total_rel_improvement:.2f}% of the previously unknown proteins.")
    
    # Statistical Significance Test
    if len(fam_df) >= 3:
        _, p_val = stats.wilcoxon(fam_df['Initial_Undesc_PCT'], fam_df['Enhanced_Undesc_PCT'])
        print(f"  - Significance: p = {p_val:.2e}")
    print("-" * 30)

# Global Summary
overall_reduction = df['Rel_Reduction_PCT'].mean()
print(f"\nOVERALL PIPELINE PERFORMANCE: {overall_reduction:.2f}% of hypothetical space resolved across all families.")



# ---------------------------------------------------------
# B. VISUALIZATION 1: Combined Hypothetical Reduction
# ---------------------------------------------------------
print("\n--- 2. Generating Combined Boxplot ---")
df_long = pd.melt(df, id_vars=['Family', 'File'], 
                  value_vars=['Initial_Undesc_PCT', 'Enhanced_Undesc_PCT'],
                  var_name='Condition', value_name='PCT')

df_long['Condition'] = df_long['Condition'].replace({
    'Initial_Undesc_PCT': 'Initial', 'Enhanced_Undesc_PCT': 'Enhanced'
})

plt.figure(figsize=(12, 7))
sns.boxplot(x='Family', y='PCT', hue='Condition', data=df_long, palette="Set2")
# Add individual points to see the distribution per genome
sns.stripplot(x='Family', y='PCT', hue='Condition', data=df_long, 
              dodge=True, color='black', alpha=0.3, legend=False)

plt.title("Reduction in Hypothetical Space across Cyanobacterial Families")
plt.ylabel("Undescribed Coding Space (%)")
plt.savefig(f"{OUTPUT_DIR}/Combined_Hypothetical_Reduction.png", dpi=300)

# ---------------------------------------------------------
# C. VISUALIZATION 2: Functional Gains (GO/COG/KEGG)
# ---------------------------------------------------------
print("\n--- 3. Generating Combined Functional Gains ---")
metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']
plot_data = []

for family in df['Family'].unique():
    fam_df = df[df['Family'] == family]
    for m in metrics:
        col_init = f'Initial_{m}'
        col_enh = f'Enhanced_{m}'
        
        if col_init in df.columns and col_enh in df.columns:
            # Average per genome
            init_val = pd.to_numeric(fam_df[col_init], errors='coerce').mean()
            enh_val = pd.to_numeric(fam_df[col_enh], errors='coerce').mean()
            
            # Use 'Metric' for X-axis and 'Family + Condition' for the bars
            label = m.split('_')[0]
            plot_data.append({'Metric': label, 'Family': family, 'Condition': 'Initial', 'Avg_Count': init_val})
            plot_data.append({'Metric': label, 'Family': family, 'Condition': 'Enhanced', 'Avg_Count': enh_val})

df_plot = pd.DataFrame(plot_data)

# 2. Create a combined Label for the legend
# This creates groups like "Nostocaceae - Initial"
df_plot['Group'] = df_plot['Family'] + " (" + df_plot['Condition'] + ")"

# 3. Plotting
plt.figure(figsize=(14, 8))
sns.barplot(
    data=df_plot, 
    x='Metric', 
    y='Avg_Count', 
    hue='Group', 
    palette='Paired' # 'Paired' works well for showing Initial/Enhanced side-by-side
)

plt.title("Functional Annotation Gains Across All Cyanobacterial Families", fontsize=16)
plt.ylabel("Average Annotations per Genome", fontsize=12)
plt.xlabel("Annotation Database", fontsize=12)
plt.legend(title="Family & Condition", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.savefig(f"{OUTPUT_DIR}/Figure_Grouped_Functional_Gains.png", dpi=300)

# ---------------------------------------------------------
# D. VISUALIZATION 3: Pseudogene Addition
# ---------------------------------------------------------
print("\n--- 4. Generating Combined Pseudogene Percentage Chart ---")

plt.figure(figsize=(12, 7))
pseudo_long = pd.melt(df, id_vars=['Family'], 
                      value_vars=['Initial_Pseudo_PCT', 'Enhanced_Pseudo_PCT'],
                      var_name='Condition', value_name='Percentage')

pseudo_long['Condition'] = pseudo_long['Condition'].replace({
    'Initial_Pseudo_PCT': 'Initial', 'Enhanced_Pseudo_PCT': 'Enhanced'
})

sns.boxplot(x='Family', y='Percentage', hue='Condition', data=pseudo_long, palette="viridis")
plt.title("Pseudogene Percentage relative to Total Estimated Genes")
plt.ylabel("Pseudogenes (%)")
plt.savefig(f"{OUTPUT_DIR}/Combined_Pseudogene_Percentage.png", dpi=300)


print(f"\nDone! Combined plots saved in '{OUTPUT_DIR}'")