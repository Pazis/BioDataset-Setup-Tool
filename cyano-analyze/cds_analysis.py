import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

# --- CONFIGURATION ---
INPUT_FILE = "/home/angelos/Desktop/Thesis/test/full-mode/goodQ_summary.csv"  # Ensure this matches your filename
OUTPUT_DIR = "plots"

if not os.path.exists(INPUT_FILE):
    print(f"Error: Could not find '{INPUT_FILE}'. Please run your Nextflow pipeline first.")
    sys.exit(1)

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# 1. Load Data
df = pd.read_csv(INPUT_FILE)

# Clean percentage strings
def clean_pct(val):
    if isinstance(val, str):
        return float(val.replace('%', ''))
    return float(val)

df['Initial_Undesc_PCT'] = df['Initial_Undesc_PCT'].apply(clean_pct)
df['Enhanced_Undesc_PCT'] = df['Enhanced_Undesc_PCT'].apply(clean_pct)

# ---------------------------------------------------------
# A. STATISTICAL ANALYSIS (Hypothetical Reduction)
# ---------------------------------------------------------
print("\n--- 1. Statistical Analysis: Hypothetical Reduction ---")

diff = df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']

if len(diff) < 3:
    print("Warning: Not enough data points for normality test (need > 3).")
    shapiro_p = 1.0 
else:
    shapiro_stat, shapiro_p = stats.shapiro(diff)

print(f"Normality Test (p-value): {shapiro_p:.5f}")

if shapiro_p > 0.05:
    print("-> Data distribution is Normal. Using Paired t-test.")
    stat, p_val = stats.ttest_rel(df['Initial_Undesc_PCT'], df['Enhanced_Undesc_PCT'])
    test_name = "Paired t-test"
else:
    print("-> Data distribution is Non-Normal. Using Wilcoxon Signed-Rank Test.")
    stat, p_val = stats.wilcoxon(df['Initial_Undesc_PCT'], df['Enhanced_Undesc_PCT'])
    test_name = "Wilcoxon Test"

print(f"Test Used: {test_name}")
print(f"Statistic: {stat:.2f}")
print(f"P-Value:   {p_val:.5e}")
print(f"Average Hypothetical Reduction: {diff.mean():.2f}%")

# ---------------------------------------------------------
# NEW: CALCULATE % INCREASE PER DATABASE
# ---------------------------------------------------------
print("\n" + "="*60)
print("FUNCTIONAL GAINS: PERCENTAGE INCREASE PER DB")
print("="*60)

metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']
pct_increases = {}

for metric in metrics:
    col_init = f'Initial_{metric}'
    col_enh = f'Enhanced_{metric}'
    
    if col_init in df.columns and col_enh in df.columns:
        # Sum counts across all samples
        total_init = pd.to_numeric(df[col_init], errors='coerce').sum()
        total_enh = pd.to_numeric(df[col_enh], errors='coerce').sum()
        
        # Calculate % Increase
        if total_init > 0:
            pct_increase = ((total_enh - total_init) / total_init) * 100
        else:
            pct_increase = 0
            
        pct_increases[metric.replace('_entries', '')] = pct_increase
        
        print(f"{metric.replace('_entries', ''):<10} | Initial: {total_init:<6} -> Enhanced: {total_enh:<6} | Increase: +{pct_increase:.2f}%")

print("="*60)

# ---------------------------------------------------------
# NEW: TOTAL IMPROVEMENT CALCULATION (Relative Reduction)
# ---------------------------------------------------------
raw_improvement = (
    (df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']) / 
    df['Initial_Undesc_PCT'].replace(0, 1)
) * 100

df['Rel_Reduction_PCT'] = raw_improvement.abs()

avg_init = df['Initial_Undesc_PCT'].mean()
avg_enh = df['Enhanced_Undesc_PCT'].mean()
avg_total_improvement = df['Rel_Reduction_PCT'].mean()

print(f"\nAverage Initial Hypotheticals:  {avg_init:.2f}%")
print(f"Average Enhanced Hypotheticals: {avg_enh:.2f}%")
print(f"TOTAL HYPOTHETICAL REDUCTION EFFICIENCY: {avg_total_improvement:.2f}%")

# ---------------------------------------------------------
# PLOTS
# ---------------------------------------------------------

# Figure 1: Paired Boxplot
print("\n--- Generating Figure 1 (Paired Boxplot) ---")
df_long = pd.melt(df, 
                  id_vars=['File'], 
                  value_vars=['Initial_Undesc_PCT', 'Enhanced_Undesc_PCT'],
                  var_name='Condition', value_name='Hypothetical_Percentage')
df_long['Condition'] = df_long['Condition'].replace({
    'Initial_Undesc_PCT': 'Initial (Bakta)', 'Enhanced_Undesc_PCT': 'Enhanced (Pipeline)'
})
plt.figure(figsize=(8, 6))
sns.violinplot(x='Condition', y='Hypothetical_Percentage', data=df_long, palette="Set2", inner=None)
sns.stripplot(x='Condition', y='Hypothetical_Percentage', data=df_long, color='black', alpha=0.5)
for i in range(len(df)):
    y1 = df['Initial_Undesc_PCT'].iloc[i]
    y2 = df['Enhanced_Undesc_PCT'].iloc[i]
    plt.plot([0, 1], [y1, y2], color='grey', linewidth=0.5, alpha=0.5)
plt.title(f"Reduction in Hypothetical Coding Space\n({test_name}, p < {p_val:.1e})")
plt.ylabel("Undescribed Coding Space (%)")
plt.savefig(f"{OUTPUT_DIR}/Figure1_Hypothetical_Reduction.png", dpi=300)

# Figure 2: Functional Gains (Counts)
print("\n--- Generating Figure 2 (Functional Gains Counts) ---")
plot_data = []
for metric in metrics:
    col_init = f'Initial_{metric}'
    col_enh = f'Enhanced_{metric}'
    if col_init in df.columns and col_enh in df.columns:
        init_sum = pd.to_numeric(df[col_init], errors='coerce').sum()
        enh_sum = pd.to_numeric(df[col_enh], errors='coerce').sum()
        plot_data.append({'Metric': metric.replace('_entries', ''), 'Condition': 'Initial', 'Count': init_sum})
        plot_data.append({'Metric': metric.replace('_entries', ''), 'Condition': 'Enhanced', 'Count': enh_sum})
if plot_data:
    df_counts = pd.DataFrame(plot_data)
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Metric', y='Count', hue='Condition', data=df_counts, palette="colorblind")
    plt.title("Total Functional Annotations Recovered")
    plt.ylabel("Total Count")
    plt.savefig(f"{OUTPUT_DIR}/Figure2_Functional_Gains.jpeg", dpi=300)

# Figure 3: Contigs vs Hypotheticals
print("\n--- Generating Figure 3 (Contigs vs Hypotheticals) ---")
if 'Contigs' in df.columns:
    df_sorted = df.sort_values(by='Contigs')
    plt.figure(figsize=(10, 6))
    plt.plot(df_sorted['Contigs'], df_sorted['Initial_Undesc_PCT'], marker='o', label='Initial', linestyle='-', alpha=0.7)
    plt.plot(df_sorted['Contigs'], df_sorted['Enhanced_Undesc_PCT'], marker='o', label='Enhanced', linestyle='-', alpha=0.7)
    plt.title("Impact of Assembly Quality on Annotations")
    plt.xlabel("Number of Contigs")
    plt.ylabel("Undescribed Coding Space (%)")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.savefig(f"{OUTPUT_DIR}/Figure3_Contigs_vs_Hypotheticals.png", dpi=300)

# Figure 4: Pseudogenes
print("\n--- Generating Figure 4 (Pseudogenes) ---")
if 'Initial_Pseudogene_candidates' in df.columns:
    pseudo_data = []
    for _, row in df.iterrows():
        pseudo_data.append({'File': row['File'], 'Condition': 'Initial', 'Count': row['Initial_Pseudogene_candidates']})
        pseudo_data.append({'File': row['File'], 'Condition': 'Enhanced', 'Count': row['Enhanced_Pseudogene_candidates']})
    df_pseudo = pd.DataFrame(pseudo_data)
    plt.figure(figsize=(12, 6))
    sns.barplot(x='File', y='Count', hue='Condition', data=df_pseudo, palette="viridis")
    plt.title("Comparison of Detected Pseudogenes")
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Number of Pseudogenes")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/Figure4_Pseudogene_Addition.png", dpi=300)

# Figure 5: Total Improvement
print("\n--- Generating Figure 5 (Total Improvement) ---")
plt.figure(figsize=(10, 6))
sns.barplot(x='File', y='Rel_Reduction_PCT', data=df, palette="viridis")
plt.axhline(avg_total_improvement, color='red', linestyle='--', label=f'Mean Improvement ({avg_total_improvement:.1f}%)')
plt.title("Pipeline Efficiency: % of Hypothetical Proteins Annotated")
plt.ylabel("Relative Reduction (%)")
plt.xlabel("Sample File")
plt.xticks(rotation=45, ha='right')
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/Figure5_Total_Improvement.png", dpi=300)

# ---------------------------------------------------------
# NEW VISUALIZATION 6: Percentage Increase per DB Bar Chart
# ---------------------------------------------------------
print("\n--- Generating Figure 6 (DB % Increase) ---")
if pct_increases:
    db_df = pd.DataFrame(list(pct_increases.items()), columns=['Database', 'Percent_Increase'])
    
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(x='Database', y='Percent_Increase', data=db_df, palette="muted")
    
    # Add text labels on top of bars
    for i in ax.containers:
        ax.bar_label(i, fmt='%.1f%%', padding=3)
        
    plt.title("Percentage Increase in Functional Annotations per Database")
    plt.ylabel("Increase (%)")
    plt.ylim(0, db_df['Percent_Increase'].max() * 1.15) # Add headroom for labels
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/Figure6_DB_Percentage_Increase.png", dpi=300)
    print("-> Saved plots/Figure6_DB_Percentage_Increase.png")

print("\nDone! Check the 'plots' folder.")