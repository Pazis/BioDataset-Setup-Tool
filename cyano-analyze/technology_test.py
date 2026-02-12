import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os

# --- CONFIGURATION ---
# Οργάνωση των αρχείων σου ανά τεχνολογία: Illumina vs PacBio
INPUT_FILES = [
    ("/home/angelos/Desktop/Thesis/test/full-mode/summary_illumina.csv", "Illumina"),
    ("/home/angelos/Desktop/Thesis/test/full-mode/summary_pacbio.csv", "PacBio")
]
OUTPUT_DIR = "technology_comparison_plots"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# --- 1. Load and Combine Data ---
all_data = []

def clean_pct(val):
    if isinstance(val, str):
        # Αφαίρεση του % και μετατροπή σε float, διατήρηση του απόλυτου αριθμού αν υπάρχει -%
        return abs(float(str(val).replace('%', '')))
    return abs(float(val))

for file_path, tech in INPUT_FILES:
    if os.path.exists(file_path):
        temp_df = pd.read_csv(file_path)
        temp_df['Technology'] = tech  # Προσθήκη στήλης για την τεχνολογία
        all_data.append(temp_df)
    else:
        print(f"Warning: Skipping {file_path}, file not found.")

if not all_data:
    print("Error: No data loaded. Check your file paths.")
    sys.exit(1)

df = pd.concat(all_data, ignore_index=True)

# Καθαρισμός δεδομένων
df['Initial_Undesc_PCT'] = df['Initial_Undesc_PCT'].apply(clean_pct)
df['Enhanced_Undesc_PCT'] = df['Enhanced_Undesc_PCT'].apply(clean_pct)

# Υπολογισμός Relative Reduction (Πόσο % των αρχικά "άγνωστων" λύθηκε)
df['Rel_Reduction_PCT'] = (
    (df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']) / 
    df['Initial_Undesc_PCT'].replace(0, 1)
) * 100
df['Rel_Reduction_PCT'] = df['Rel_Reduction_PCT'].abs() # Μετατροπή τυχόν αρνητικών σε θετικά

# ---------------------------------------------------------
# A. STATISTICAL ANALYSIS (Per Technology)
# ---------------------------------------------------------
print("\n" + "="*60)
print("STATISTICAL ANALYSIS BY SEQUENCING TECHNOLOGY")
print("="*60)

for tech in df['Technology'].unique():
    tech_df = df[df['Technology'] == tech]
    
    avg_initial = tech_df['Initial_Undesc_PCT'].mean()
    avg_enhanced = tech_df['Enhanced_Undesc_PCT'].mean()
    total_rel_improvement = tech_df['Rel_Reduction_PCT'].mean()
    
    print(f"TECHNOLOGY: {tech}")
    print(f"  - Avg Initial Hypothetical:  {avg_initial:.2f}%")
    print(f"  - Avg Enhanced Hypothetical: {avg_enhanced:.2f}%")
    print(f"  - PIPELINE IMPACT: Resolved {total_rel_improvement:.2f}% of unknowns.")
    
    if len(tech_df) >= 3:
        _, p_val = stats.wilcoxon(tech_df['Initial_Undesc_PCT'], tech_df['Enhanced_Undesc_PCT'])
        print(f"  - Statistical Significance: p = {p_val:.2e}")
    print("-" * 30)

# ---------------------------------------------------------
# B. VISUALIZATION 1: Boxplot ανά Τεχνολογία
# ---------------------------------------------------------
print("\n--- 2. Generating Comparison Boxplot ---")
df_long = pd.melt(df, id_vars=['Technology', 'File'], 
                  value_vars=['Initial_Undesc_PCT', 'Enhanced_Undesc_PCT'],
                  var_name='Condition', value_name='PCT')

df_long['Condition'] = df_long['Condition'].replace({
    'Initial_Undesc_PCT': 'Initial (Bakta)', 'Enhanced_Undesc_PCT': 'Enhanced (Pipeline)'
})

plt.figure(figsize=(10, 6))
sns.boxplot(x='Technology', y='PCT', hue='Condition', data=df_long, palette="Set2")
sns.stripplot(x='Technology', y='PCT', hue='Condition', data=df_long, 
              dodge=True, color='black', alpha=0.3, legend=False)

plt.title("Hypothetical Space Reduction: Illumina vs PacBio")
plt.ylabel("Undescribed Coding Space (%)")
plt.savefig(f"{OUTPUT_DIR}/Technology_Reduction_Comparison.png", dpi=300)

# ---------------------------------------------------------
# NEW: CALCULATE % INCREASE PER DB BY TECHNOLOGY
# ---------------------------------------------------------
print("\n" + "="*60)
print("FUNCTIONAL GAINS: PERCENTAGE INCREASE PER DB (BY TECH)")
print("="*60)

metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']
db_increase_data = []

for tech in df['Technology'].unique():
    tech_df = df[df['Technology'] == tech]
    print(f"TECHNOLOGY: {tech}")
    
    for metric in metrics:
        col_init = f'Initial_{metric}'
        col_enh = f'Enhanced_{metric}'
        
        if col_init in df.columns and col_enh in df.columns:
            # Sum counts across all samples in that technology group
            total_init = pd.to_numeric(tech_df[col_init], errors='coerce').sum()
            total_enh = pd.to_numeric(tech_df[col_enh], errors='coerce').sum()
            
            # Calculate % Increase
            if total_init > 0:
                pct_increase = ((total_enh - total_init) / total_init) * 100
            else:
                pct_increase = 0
            
            clean_metric = metric.replace('_entries', '')
            print(f"  {clean_metric:<5} | Init: {total_init:<6} -> Enh: {total_enh:<6} | Increase: +{pct_increase:.2f}%")
            
            # Save for plotting
            db_increase_data.append({
                'Technology': tech,
                'Database': clean_metric,
                'Percent_Increase': pct_increase
            })
    print("-" * 30)
print("="*60)

# ---------------------------------------------------------
# B. VISUALIZATION 1: Boxplot ανά Τεχνολογία
# ---------------------------------------------------------
print("\n--- Generating Figure 1 (Comparison Boxplot) ---")
df_long = pd.melt(df, id_vars=['Technology', 'File'], 
                  value_vars=['Initial_Undesc_PCT', 'Enhanced_Undesc_PCT'],
                  var_name='Condition', value_name='PCT')

df_long['Condition'] = df_long['Condition'].replace({
    'Initial_Undesc_PCT': 'Initial', 'Enhanced_Undesc_PCT': 'Enhanced'
})

plt.figure(figsize=(10, 6))
sns.boxplot(x='Technology', y='PCT', hue='Condition', data=df_long, palette="Set2")
sns.stripplot(x='Technology', y='PCT', hue='Condition', data=df_long, 
              dodge=True, color='black', alpha=0.3, legend=False)

plt.title("Hypothetical Space Reduction: Illumina vs PacBio")
plt.ylabel("Undescribed Coding Space (%)")
plt.savefig(f"{OUTPUT_DIR}/Figure1_Technology_Reduction.png", dpi=300)

# ---------------------------------------------------------
# C. VISUALIZATION 2: Functional Gains per Technology
# ---------------------------------------------------------
print("\n--- 3. Generating Functional Gains Chart ---")
metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']
plot_data = []

for tech in df['Technology'].unique():
    tech_df = df[df['Technology'] == tech]
    for m in metrics:
        col_init = f'Initial_{m}'
        col_enh = f'Enhanced_{m}'
        if col_init in df.columns and col_enh in df.columns:
            init_val = pd.to_numeric(tech_df[col_init], errors='coerce').mean()
            enh_val = pd.to_numeric(tech_df[col_enh], errors='coerce').mean()
            
            label = m.split('_')[0]
            plot_data.append({'Metric': label, 'Technology': tech, 'Condition': 'Initial', 'Avg_Count': init_val})
            plot_data.append({'Metric': label, 'Technology': tech, 'Condition': 'Enhanced', 'Avg_Count': enh_val})

df_plot = pd.DataFrame(plot_data)
df_plot['Group'] = df_plot['Technology'] + " - " + df_plot['Condition']

plt.figure(figsize=(12, 7))
sns.barplot(data=df_plot, x='Metric', y='Avg_Count', hue='Group', palette='Paired')
plt.title("Average Functional Annotation Recovery by Technology")
plt.ylabel("Mean Annotations per Genome")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/Technology_Functional_Gains.png", dpi=300)

# ---------------------------------------------------------
# D. VISUALIZATION 3: Impact Magnitude
# ---------------------------------------------------------
print("\n--- 4. Generating Efficiency Comparison ---")
plt.figure(figsize=(8, 6))
sns.barplot(x='Technology', y='Rel_Reduction_PCT', data=df, palette="muted", errorbar=None)
plt.title("Pipeline Efficiency: % of Hypotheticals Resolved")
plt.ylabel("Relative Improvement (%)")
plt.savefig(f"{OUTPUT_DIR}/Technology_Efficiency_Comparison.png", dpi=300)

# ---------------------------------------------------------
# E. VISUALIZATION 4: Percentage Increase by Technology (NEW)
# ---------------------------------------------------------
print("\n--- Generating Figure 4 (Percentage Increase by Tech) ---")
if db_increase_data:
    df_db_inc = pd.DataFrame(db_increase_data)
    
    plt.figure(figsize=(10, 6))
    
    # Plot: X=Database, Y=Increase%, Hue=Technology
    sns.barplot(data=df_db_inc, x='Database', y='Percent_Increase', hue='Technology', palette='muted')
    
    # Add labels on top of bars
    ax = plt.gca()
    for i in ax.containers:
        ax.bar_label(i, fmt='%.1f%%', padding=3, fontsize=10)
        
    plt.title("Percentage Increase in Functional Annotations by Technology", fontsize=14)
    plt.ylabel("Increase in Annotations (%)")
    plt.xlabel("Annotation Database")
    plt.legend(title="Sequencing Tech")
    plt.tight_layout()
    
    plt.savefig(f"{OUTPUT_DIR}/Figure4_DB_Percentage_Increase_By_Tech.png", dpi=300)



print(f"\nDone! Plots saved in '{OUTPUT_DIR}'")

