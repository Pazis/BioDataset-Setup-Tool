import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import glob 

# ---------------------------------------------------------
# 1. DATA LOADING & PREP
# ---------------------------------------------------------
def count_unique_contigs(gff_path):
    """Counts unique sequence regions (contigs) in a GFF file."""
    unique_contigs = set()
    try:
        with open(gff_path, 'r') as file:
            for line in file:
                if line.startswith('#') or not line.strip():
                    continue
                columns = line.split('\t')
                if len(columns) > 0:
                    unique_contigs.add(columns[0])
        return len(unique_contigs)
    except Exception:
        return 0

def calculate_genome_size(gff_path):
    """Calculates total genome size (bp) from a GFF file."""
    contig_sizes = {}
    try:
        with open(gff_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line: continue
                
                # Method 1: Header (Accurate)
                if line.startswith('##sequence-region'):
                    parts = line.split()
                    if len(parts) >= 4:
                        contig_sizes[parts[1]] = int(parts[3])
                    continue

                if line.startswith('#'): continue
                
                # Method 2: Features (Fallback)
                cols = line.split('\t')
                if len(cols) >= 5:
                    seqid = cols[0]
                    try:
                        end = int(cols[4])
                        # Keep the largest coordinate seen for this contig
                        if end > contig_sizes.get(seqid, 0):
                            contig_sizes[seqid] = end
                    except ValueError:
                        continue
        return sum(contig_sizes.values())
    except Exception:
        return 0

def load_csv(file_path):
    df = pd.read_csv(file_path)
    def clean_pct(val):
        if isinstance(val, str):
            return float(val.replace('%', ''))
        return float(val)

    # Clean percentage columns for all 3 tools
    if 'Initial_Undesc_PCT' in df.columns:
        df['Initial_Undesc_PCT'] = df['Initial_Undesc_PCT'].apply(clean_pct)
    if 'Enhanced_Undesc_PCT' in df.columns:
        df['Enhanced_Undesc_PCT'] = df['Enhanced_Undesc_PCT'].apply(clean_pct)
    if 'Prokka_Undesc_PCT' in df.columns:
        df['Prokka_Undesc_PCT'] = df['Prokka_Undesc_PCT'].apply(clean_pct)
        
    return df

def enrich_data(df, gff_dir):
    """Adds 'Contigs' and 'Genome_Size' columns from GFF files."""
    print("\n--- DEBUG: Calculating Contigs & Genome Size ---")
    print(f"1. Looking in directory: '{gff_dir}'")
    
    if not gff_dir or not os.path.isdir(gff_dir):
        print("   [ERROR] Directory not found.")
        return df

    # Map files
    all_gffs = glob.glob(os.path.join(gff_dir, "*"))
    gff_files = [f for f in all_gffs if f.endswith('.gff') or f.endswith('.gff3')]
    print(f"2. Found {len(gff_files)} GFF files.")
    
    gff_map = {}
    for p in gff_files:
        base = os.path.basename(p)
        # Strategy A: Exact Name
        gff_map[base] = p
        # Strategy B: No extension
        gff_map[os.path.splitext(base)[0]] = p
        
        # Strategy C: Clean suffix (handle _enhanced)
        clean = base.replace(".gff3", "").replace(".gff", "")
        if clean.endswith("_enhanced"): 
            clean = clean.replace("_enhanced", "")
        gff_map[clean] = p

    contigs = []
    sizes = []
    matches = 0
    
    # Iterate through CSV rows
    for index, row in df.iterrows():
        csv_name = str(row['File']).strip()
        target = gff_map.get(csv_name)
        
        if target:
            contigs.append(count_unique_contigs(target))
            sizes.append(calculate_genome_size(target))
            matches += 1
        else:
            if matches == 0 and index < 3:
                print(f"   [FAIL] CSV says '{csv_name}' -> No match in GFF map.")
            contigs.append(0)
            sizes.append(0)

    df['Contigs'] = contigs
    df['Genome_Size'] = sizes
    # Avoid division by zero for empty files
    df['Genome_Size_MB'] = df['Genome_Size'].apply(lambda x: x / 1_000_000 if x > 0 else 0)
    
    print(f"4. Result: Enriched {matches}/{len(df)} files.")
    return df

# ---------------------------------------------------------
# 2. STATISTICAL FUNCTIONS (UNCHANGED - Initial vs Enhanced Only)
# ---------------------------------------------------------
def correlate_size(df):
    """
    Correlates Genome Size with Hypothetical % using Pearson (if normal) 
    or Spearman (if non-normal).
    """
    print("\n--- 2. Genome Size Correlation ---")
    if 'Genome_Size' not in df.columns or df['Genome_Size'].sum() == 0:
        print("Skipping: No size data.")
        return

    # 1. Normality Checks (Shapiro-Wilk)
    stat_size, p_size = stats.shapiro(df['Genome_Size'])
    stat_hypo, p_hypo = stats.shapiro(df['Enhanced_Undesc_PCT'])

    print(f"   Normality (Genome Size): p={p_size:.5f}")
    print(f"   Normality (Hypothetical %): p={p_hypo:.5f}")

    # 2. Choose Test
    if p_size > 0.05 and p_hypo > 0.05:
        print("-> Both Normal. Using Pearson Correlation.")
        test_type = "Pearson (r)"
        corr, p_val = stats.pearsonr(df['Genome_Size'], df['Enhanced_Undesc_PCT'])
    else:
        print("-> Non-Normal Data. Using Spearman Correlation.")
        test_type = "Spearman (rho)"
        corr, p_val = stats.spearmanr(df['Genome_Size'], df['Enhanced_Undesc_PCT'])
    
    # 3. Report Results
    print(f"   Correlation: {test_type}={corr:.3f}, p={p_val:.5e}")
    
    if p_val < 0.05:
        strength = "Strong" if abs(corr) > 0.66 else "Moderate" if abs(corr) > 0.33 else "Weak"
        direction = "Positive" if corr > 0 else "Negative"
        print(f"-> Significant {strength} {direction} correlation detected.")
    else:
        print("-> No significant correlation.")

def hypothetical_reduction(df):
    """Performs normality test and paired t-test/Wilcoxon on reduction data (Initial vs Enhanced)."""
    print("\n--- 1. Statistical Analysis: Hypothetical Reduction (Initial vs Enhanced) ---")

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
    
    return test_name, p_val, diff.mean()

def per_database_increase(df):
    """Calculates percentage increase for specific database columns (Initial vs Enhanced)."""
    print("\n" + "="*60)
    print("FUNCTIONAL GAINS: PERCENTAGE INCREASE PER DB (Initial -> Enhanced)")
    print("="*60)

    metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']
    pct_increases = {}

    for metric in metrics:
        col_init = f'Initial_{metric}'
        col_enh = f'Enhanced_{metric}'
        
        if col_init in df.columns and col_enh in df.columns:
            total_init = pd.to_numeric(df[col_init], errors='coerce').sum()
            total_enh = pd.to_numeric(df[col_enh], errors='coerce').sum()
            
            if total_init > 0:
                pct_increase = ((total_enh - total_init) / total_init) * 100
            else:
                pct_increase = 0
                
            pct_increases[metric.replace('_entries', '')] = pct_increase
            
            print(f"{metric.replace('_entries', ''):<10} | Initial: {total_init:<6} -> Enhanced: {total_enh:<6} | Increase: +{pct_increase:.2f}%")

    print("="*60)
    return pct_increases

def total_improvement(df):
    """Calculates the relative reduction percentage (Initial vs Enhanced)."""
    # Avoid division by zero
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
    
    return avg_total_improvement

# ---------------------------------------------------------
# 3. PLOTTING FUNCTION (MODIFIED TO INCLUDE PROKKA)
# ---------------------------------------------------------

def generate_plots(df, test_name, p_val, pct_increases, avg_total_improvement, OUTPUT_DIR):
    """Generates Figures 1 through 7, including Prokka in visualizations."""
    
    metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']

    TOOL_COLORS = {
    'Prokka': '#0072B2',   # Strong Blue
    'Initial': '#999999',  # Grey
    'Enhanced': '#D55E00'  # Vermilion
}
    

    # -------------------------------------------------------
    # Figure 1: Comparison Boxplot (Prokka vs Initial vs Enhanced)
    # -------------------------------------------------------
    print("\n--- Generating Figure 1 (Comparison Boxplot) ---")
    
    # We melt all 3 columns now
    melt_vars = ['Initial_Undesc_PCT', 'Enhanced_Undesc_PCT']
    if 'Prokka_Undesc_PCT' in df.columns:
        melt_vars.insert(0, 'Prokka_Undesc_PCT') # Add Prokka first
        
    df_long = pd.melt(df, 
                      id_vars=['File'], 
                      value_vars=melt_vars,
                      var_name='Condition', value_name='Hypothetical_Percentage')
    
    df_long['Condition'] = df_long['Condition'].replace({
        'Prokka_Undesc_PCT': 'Prokka',
        'Initial_Undesc_PCT': 'Initial (Bakta)', 
        'Enhanced_Undesc_PCT': 'Enhanced (Pipeline)'
    })
    
    plt.figure(figsize=(9, 6))
    # Using a 3-color palette
    sns.boxplot(x='Condition', y='Hypothetical_Percentage', data=df_long, palette=[TOOL_COLORS['Prokka'], TOOL_COLORS['Initial'], TOOL_COLORS['Enhanced']])
    
    # Title notes the significance between Init/Enhanced, but plot shows all 3
    plt.title(f"Comparison of Undescribed Coding Space\n(Init vs Enhanced: {test_name}, p < {p_val:.1e})")
    plt.ylabel("Undescribed Coding Space (%)")
    plt.savefig(f"{OUTPUT_DIR}/Figure1_Hypothetical_Reduction.png", dpi=300)

    # -------------------------------------------------------
    # Figure 2: Functional Gains (Counts)
    # -------------------------------------------------------
    print("\n--- Generating Figure 2 (Functional Gains Counts) ---")
    plot_data = []
    
    # List of prefixes to check
    prefixes = [('Prokka', 'Prokka'), ('Initial', 'Initial'), ('Enhanced', 'Enhanced')]
    
    for metric in metrics:
        clean_metric_name = metric.replace('_entries', '')
        
        for prefix, label in prefixes:
            col_name = f'{prefix}_{metric}'
            if col_name in df.columns:
                total_sum = pd.to_numeric(df[col_name], errors='coerce').sum()
                plot_data.append({'Metric': clean_metric_name, 'Condition': label, 'Count': total_sum})
                
    if plot_data:
        df_counts = pd.DataFrame(plot_data)
        plt.figure(figsize=(12, 6))
        sns.barplot(x='Metric', y='Count', hue='Condition', data=df_counts, palette=[TOOL_COLORS['Prokka'],TOOL_COLORS['Initial'],TOOL_COLORS['Enhanced']])
        plt.title("Total Functional Annotations Recovered")
        plt.ylabel("Total Count")
        plt.legend(title="Tool")
        plt.savefig(f"{OUTPUT_DIR}/Figure2_Functional_Gains.jpeg", dpi=300)

    # -------------------------------------------------------
    # Figure 3: Contigs vs Hypotheticals (Scatter with Prokka)
    # -------------------------------------------------------
    print("\n--- Generating Figure 3 (Contigs vs Hypotheticals) ---")
    
    if 'Contigs' in df.columns and df['Contigs'].sum() > 0:
        df_sorted = df.sort_values(by='Contigs')
        
        plt.figure(figsize=(10, 6))
        
        # 1. Prokka (Blue) - if exists
        if 'Prokka_Undesc_PCT' in df.columns:
            plt.plot(df_sorted['Contigs'], df_sorted['Prokka_Undesc_PCT'], 
                     marker='^', label='Prokka', linestyle='--', color=TOOL_COLORS['Prokka'], alpha=0.6)

        # 2. Initial (Grey)
        plt.plot(df_sorted['Contigs'], df_sorted['Initial_Undesc_PCT'], 
                 marker='o', label='Initial (Bakta)', linestyle='-', color=TOOL_COLORS['Initial'], alpha=0.7)
        
        # 3. Enhanced (Red)
        plt.plot(df_sorted['Contigs'], df_sorted['Enhanced_Undesc_PCT'], 
                 marker='o', label='Enhanced', linestyle='-', color=TOOL_COLORS['Enhanced'], alpha=0.8)
        
        plt.title("Impact of Assembly Quality on Annotations")
        plt.xlabel("Number of Contigs")
        plt.ylabel("Undescribed Coding Space (%)")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(f"{OUTPUT_DIR}/Figure3_Contigs_vs_Hypotheticals.png", dpi=300)
    else:
        print("Skipping Figure 3: No GFF directory provided or no contigs found.")

    # -------------------------------------------------------
    # Figure 4: Pseudogenes (Bar Plot)
    # -------------------------------------------------------
    print("\n--- Generating Figure 4 (Pseudogenes) ---")
    
    pseudo_cols = [
        ('Prokka_Pseudogene_candidates', 'Prokka'),
        ('Initial_Pseudogene_candidates', 'Initial'), 
        ('Enhanced_Pseudogene_candidates', 'Enhanced')
    ]
    
    pseudo_data = []
    for _, row in df.iterrows():
        for col, label in pseudo_cols:
            if col in df.columns:
                pseudo_data.append({'File': row['File'], 'Condition': label, 'Count': row[col]})

    if pseudo_data:
        df_pseudo = pd.DataFrame(pseudo_data)
        plt.figure(figsize=(12, 6))
        # Use simple barplot or grouped barplot
        sns.barplot(x='File', y='Count', hue='Condition', data=df_pseudo, palette="viridis")
        plt.title("Comparison of Detected Pseudogenes")
        plt.xticks([]) # Hide x-axis labels if too many genomes
        plt.ylabel("Number of Pseudogenes")
        plt.tight_layout()
        plt.savefig(f"{OUTPUT_DIR}/Figure4_Pseudogene_Addition.png", dpi=300)

    # -------------------------------------------------------
    # Figure 5: Total Improvement (Initial vs Enhanced ONLY)
    # -------------------------------------------------------
    print("\n--- Generating Figure 5 (Total Improvement - Init vs Enh) ---")
    # This plot remains strictly about the pipeline efficiency (Bakta vs Enhanced)
    if 'Rel_Reduction_PCT' not in df.columns:
         raw_improvement = ((df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']) / df['Initial_Undesc_PCT'].replace(0, 1)) * 100
         df['Rel_Reduction_PCT'] = raw_improvement.abs()

    plt.figure(figsize=(10, 6))
    sns.barplot(x='File', y='Rel_Reduction_PCT', data=df, palette="viridis")
    plt.axhline(avg_total_improvement, color='red', linestyle='--', label=f'Mean Improvement ({avg_total_improvement:.1f}%)')
    plt.title("Pipeline Efficiency: % of Hypothetical Proteins Annotated\n(Initial -> Enhanced)")
    plt.ylabel("Relative Reduction (%)")
    plt.xlabel("Genomes")
    plt.xticks([]) 
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/Figure5_Total_Improvement.png", dpi=300)

    # -------------------------------------------------------
    # Figure 6: DB Percentage Increase (Initial vs Enhanced ONLY)
    # -------------------------------------------------------
    print("\n--- Generating Figure 6 (DB % Increase - Init vs Enh) ---")
    if pct_increases:
        db_df = pd.DataFrame(list(pct_increases.items()), columns=['Database', 'Percent_Increase'])
        plt.figure(figsize=(8, 6))
        ax = sns.barplot(x='Database', y='Percent_Increase', data=db_df, palette="colorblind")
        for i in ax.containers:
            ax.bar_label(i, fmt='%.1f%%', padding=3)
        plt.title("Percentage Increase in Annotations (Initial -> Enhanced)")
        plt.ylabel("Increase (%)")
        plt.ylim(0, db_df['Percent_Increase'].max() * 1.15) 
        plt.tight_layout()
        plt.savefig(f"{OUTPUT_DIR}/Figure6_DB_Percentage_Increase.png", dpi=300)

    print(f"\nDone! Check the '{OUTPUT_DIR}' folder.")

    # -------------------------------------------------------
    # Figure 7: Genome Size Correlation (Enhanced only)
    # -------------------------------------------------------
    if 'Genome_Size_MB' in df.columns and df['Genome_Size_MB'].sum() > 0:
        print("\n--- Generating Figure 7 (Genome Size Correlation) ---")
        plt.figure(figsize=(10, 6))
        sns.regplotplot(x='Genome_Size_MB', y='Enhanced_Undesc_PCT', data=df, 
                    scatter_kws={'color': '#D55E00', 'label': 'Enhanced'}, 
                    line_kws={'color': 'darkred'})
        
        plt.title("Genome Size vs. Remaining Hypothetical Proteins")
        plt.xlabel("Genome Size (Mb)")
        plt.ylabel("Enhanced Undescribed Coding Space (%)")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.3)
        plt.savefig(f"{OUTPUT_DIR}/Figure7_Size_vs_Hypotheticals.png", dpi=300)
        plt.close()

# ---------------------------------------------------------
# 4. MODULE ENTRY POINT (To be called from main.py)
# ---------------------------------------------------------

def cmd_plot(args):
    """
    Main entry point for this module.
    """
    print("\n--- DEBUG: Starting Plot Command ---")
    print(f"1. Arguments received: {args}")
    
    input_file = args.input
    output_dir = args.output

    # Check if gff_dir exists in args
    gff_dir_arg = getattr(args, 'gff_dir', None)
    print(f"2. GFF Directory Argument is: '{gff_dir_arg}'")

    if not os.path.exists(input_file):
        print(f"Error: File {input_file} not found.")
        return
    
    df = load_csv(input_file)

    # STRICT CHECK: If the argument exists, run the enrichment
    if gff_dir_arg:
        print("3. Flag detected. Running enrichment...")
        df = enrich_data(df, gff_dir_arg)
    else:
        print("3. [WARNING] No GFF directory detected (-g). Skipping contig counting.")

    # 2. Run Analysis (Strictly Initial vs Enhanced)
    test_name, p_val, _ = hypothetical_reduction(df)
    pct_increases = per_database_increase(df)
    avg_total_improvement = total_improvement(df)
    size_correlation = correlate_size(df)

    # 3. Create Output Directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 4. Generate Plots (Includes Prokka)
    generate_plots(df, test_name, p_val, pct_increases, avg_total_improvement, output_dir)