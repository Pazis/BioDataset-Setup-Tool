import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# ==========================================
# CONFIGURATION
# ==========================================
# Enter the name of your Excel file here
INPUT_FILE = "/home/angelos/Desktop/Thesis/test/full-mode/bgc_family_analysis/analysis_all.xlsx"
# ==========================================

def analyze_bgc_data(file_path):
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        return

    # 1. Load Data
    df = pd.read_excel(file_path)

    # Filter out the "AVERAGE" row if it exists (so it doesn't skew statistics)
    df_clean = df[df['Locus Tag'] != 'AVERAGE'].copy()
    
    print(f"Loaded {len(df_clean)} genomes for analysis.\n")

    initial_count = len(df_clean)
    
    # Keep only rows where:
    # 1. AntiSMASH found at least one cluster
    # 2. GECCO found at least one cluster
    # 3. There is at least ONE overlapping cluster between them
    df_clean = df_clean[
        (df_clean['AntiSMASH Count'] > 0) & 
        (df_clean['GECCO Count'] > 0) &
        (df_clean['Overlapping Clusters'] > 0) 
    ].copy()
    
    dropped_count = initial_count - len(df_clean)
    if dropped_count > 0:
        print(f"-> Dropped {dropped_count} genomes because they had 0 counts or 0 overlap.")
    # ---------------------------------------------------------

    # ---------------------------------------------------------
    # PART A: Statistical Tests
    # ---------------------------------------------------------
    print("--- STATISTICAL ANALYSIS ---")

    # 1. Descriptive Statistics
    stats_summary = df_clean[['AntiSMASH Count', 'GECCO Count', 'Consolidation Ratio']].describe()
    print("\n1. Descriptive Statistics:")
    print(stats_summary)

    # 2. Normality Check (Shapiro-Wilk) to decide on test type
    # If p < 0.05, data is NOT normal.
    shapiro_as = stats.shapiro(df_clean['AntiSMASH Count'])
    shapiro_gc = stats.shapiro(df_clean['GECCO Count'])
    
    # 3. Pairwise Comparison (Sensitivity)
    # We use Wilcoxon Signed-Rank Test (non-parametric) because N is usually small (<30)
    # H0: There is no difference in the number of BGCs detected.
    w_stat, p_value = stats.wilcoxon(df_clean['AntiSMASH Count'], df_clean['GECCO Count'])
    
    print("\n2. Detection Power Comparison (Wilcoxon Signed-Rank Test):")
    print(f"   Statistic: {w_stat}")
    print(f"   P-value: {p_value:.5e}")
    if p_value < 0.05:
        print("   Result: Significant difference (AntiSMASH detects significantly more clusters).")
    else:
        print("   Result: No significant difference.")

    # 4. Correlation Analysis (Do they find BGCs in the same genomes?)
    # Spearman Correlation (Rank-based)
    corr, corr_p = stats.spearmanr(df_clean['AntiSMASH Count'], df_clean['GECCO Count'])
    print("\n3. Tool Correlation (Spearman's Rho):")
    print(f"   Correlation Coefficient: {corr:.3f}")
    print(f"   P-value: {corr_p:.5f}")
    if corr > 0.6:
        print("   Result: Strong positive correlation (Genome richness is consistently reflected).")
    else:
        print("   Result: Weak or no correlation.")

    # ---------------------------------------------------------
    # PART B: Visualization (Thesis Quality Plots)
    # ---------------------------------------------------------
    print("\n--- GENERATING PLOTS ---")
    
    # Set global style
    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 12})

    # Plot 1: Detection Sensitivity (Box Plot)
    plt.figure(figsize=(8, 6))
    
    # Prepare data for plotting
    plot_data = df_clean.melt(value_vars=['AntiSMASH Count', 'GECCO Count'], 
                              var_name='Tool', value_name='Count')
    
    sns.boxplot(x='Tool', y='Count', data=plot_data, palette=['#3498db', '#e74c3c'], width=0.5)
    sns.stripplot(x='Tool', y='Count', data=plot_data, color='black', alpha=0.3, jitter=True) # Add individual points
    
    plt.title('BGC Detection Sensitivity:\nAntiSMASH vs GECCO', fontsize=14, fontweight='bold')
    plt.ylabel('Number of BGCs per Genome')
    plt.xlabel('')
    plt.tight_layout()
    plt.savefig('Figure_1_Detection_Sensitivity.png', dpi=300)
    print("-> Saved 'Figure_1_Detection_Sensitivity.png'")

    # Plot 2: Fragmentation Analysis (Consolidation Ratio Histogram)
    plt.figure(figsize=(8, 6))
    sns.histplot(df_clean['Consolidation Ratio'], bins=10, kde=True, color='#9b59b6', edgecolor='black')
    
    # Add mean line
    mean_ratio = df_clean['Consolidation Ratio'].mean()
    plt.axvline(mean_ratio, color='red', linestyle='--', label=f'Mean: {mean_ratio:.2f}')
    
    plt.title('Locus Fragmentation Analysis\n(Consolidation Ratio Distribution)', fontsize=14, fontweight='bold')
    plt.xlabel('Ratio (AntiSMASH Clusters / GECCO Clusters)')
    plt.ylabel('Frequency (Genomes)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('Figure_2_Consolidation_Ratio.png', dpi=300)
    print("-> Saved 'Figure_2_Consolidation_Ratio.png'")

    # Plot 3: Boundary Precision (Scatter Plot)
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='% AntiSMASH Covered', y='% GECCO Covered', data=df_clean, 
                    s=100, color='#2ecc71', edgecolor='black', alpha=0.8)
    
    # Add diagonal line (x=y)
    plt.plot([0, 100], [0, 100], 'r--', alpha=0.5, label='Perfect Match (1:1)')
    
    plt.title('Boundary Definition Precision', fontsize=14, fontweight='bold')
    plt.xlabel('% AntiSMASH Locus Covered by GECCO (Selectivity)')
    plt.ylabel('% GECCO Locus Covered by AntiSMASH (Accuracy)')
    plt.xlim(0, 105)
    plt.ylim(50, 105) # Adjusted limit as GECCO is usually high
    plt.legend()
    plt.tight_layout()
    plt.savefig('Figure_3_Boundary_Precision.png', dpi=300)
    print("-> Saved 'Figure_3_Boundary_Precision.png'")

    plt.figure(figsize=(10, 6))
    
    # Calculate means
    avg_consensus = df_clean['Inter-tool Consensus'].mean()
    avg_efficiency = df_clean['Avg Overlap Efficiency'].mean()
    avg_ratio = df_clean['Consolidation Ratio'].mean()
    
    # Create normalized ratio for visualization (Ratio - 1, so 1.0 becomes 0 baseline)
    # Or just plot raw values. Let's plot raw values but scale Ratio to be comparable or plot on secondary axis.
    # Actually, let's keep it simple: Bar chart of averages.
    
    metrics = ['Inter-tool Consensus', 'Overlap Efficiency', 'Consolidation Ratio']
    values = [avg_consensus, avg_efficiency, avg_ratio]
    colors = ['#2ecc71', '#f1c40f', '#9b59b6'] # Green, Yellow, Purple
    
    bars = plt.bar(metrics, values, color=colors, edgecolor='black', alpha=0.8)
    
    # Add value labels on top
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height:.2f}', ha='center', va='bottom', fontweight='bold')

    plt.title('Average Pipeline Performance Metrics', fontsize=14, fontweight='bold')
    plt.ylabel('Metric Value')
    plt.ylim(0, max(values) * 1.2) # Add some headroom
    
    # Add an explanatory note on the plot
    plt.figtext(0.5, 0.01, 
                "Note: Consensus & Efficiency range [0-1]. Ratio > 1 indicates fragmentation.", 
                ha="center", fontsize=10, style='italic')

    plt.tight_layout()
    plt.savefig('Figure_4_Performance_Metrics.png', dpi=300)
    print("-> Saved 'Figure_4_Performance_Metrics.png'")

    print("\nAnalysis Complete.")

if __name__ == "__main__":
    analyze_bgc_data(INPUT_FILE)