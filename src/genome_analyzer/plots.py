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
    except Exception as e:
        print(f"Warning: Could not count contigs for {gff_path}: {e}")
        return 0
    
def load_csv(file_path):
    """Loads the CSV and cleans percentage columns."""
    df = pd.read_csv(file_path)

    def clean_pct(val):
        if isinstance(val, str):
            return float(val.replace('%', ''))
        return float(val)

    if 'Initial_Undesc_PCT' in df.columns:
        df['Initial_Undesc_PCT'] = df['Initial_Undesc_PCT'].apply(clean_pct)
    if 'Enhanced_Undesc_PCT' in df.columns:
        df['Enhanced_Undesc_PCT'] = df['Enhanced_Undesc_PCT'].apply(clean_pct)
    return df

def enrich_with_contigs(df, gff_dir):
    print(f"\n--- DEBUG: Calculating Contig Counts ---")
    print(f"1. Looking in directory: '{gff_dir}'")
    
    if not gff_dir or not os.path.isdir(gff_dir):
        print(f"   [ERROR] Directory not found or empty argument.")
        return df

    # 1. Map all available GFFs
    all_gffs = glob.glob(os.path.join(gff_dir, "*")) # Get EVERYTHING to see what's there
    gff_files = [f for f in all_gffs if f.endswith('.gff') or f.endswith('.gff3')]
    
    print(f"2. Found {len(gff_files)} GFF files.")
    if len(gff_files) > 0:
        print(f"   Example file from folder: '{os.path.basename(gff_files[0])}'")
    else:
        print("   [ERROR] No .gff or .gff3 files found! Check extensions.")
        return df

    # Create Mapping Dictionary
    gff_map = {}
    for p in gff_files:
        base = os.path.basename(p)
        
        # Strategy A: Exact Name
        gff_map[base] = p
        
        # Strategy B: Remove extensions
        no_ext = os.path.splitext(base)[0] # removes .gff3
        gff_map[no_ext] = p
        
        # Strategy C: Handle "_enhanced" suffix
        # matches: "genome_enhanced.gff3" -> "genome"
        clean_name = base.replace(".gff3", "").replace(".gff", "")
        if clean_name.endswith("_enhanced"):
            clean_name = clean_name.replace("_enhanced", "")
        gff_map[clean_name] = p

    # 2. Match CSV to GFF
    print(f"3. Matching CSV entries ({len(df)} rows) to GFF files...")
    if len(df) > 0:
        print(f"   Example entry from CSV: '{df.iloc[0]['File']}'")

    contig_counts = []
    matches = 0
    
    for index, row in df.iterrows():
        csv_name = str(row['File']).strip()
        
        # Try finding a match
        target_path = gff_map.get(csv_name)
        
        # Try extensions
        if not target_path:
            target_path = gff_map.get(f"{csv_name}.gff")
        if not target_path:
            target_path = gff_map.get(f"{csv_name}.gff3")
        # Try lowercase
        if not target_path:
            target_path = gff_map.get(csv_name.lower())

        if target_path:
            count = count_unique_contigs(target_path)
            contig_counts.append(count)
            matches += 1
        else:
            if index < 3: # Print first 3 failures
                print(f"   [FAIL] CSV says '{csv_name}', but no match found in GFF folder.")
            contig_counts.append(0)

    df['Contigs'] = contig_counts
    print(f"4. Result: Matched {matches}/{len(df)} files.")
    
    if matches == 0:
        print("   [WARNING] Total failure to match. Your CSV names likely don't match GFF filenames.")
        
    return df

# ---------------------------------------------------------
# 2. STATISTICAL FUNCTIONS
# ---------------------------------------------------------

def hypothetical_reduction(df):
    """Performs normality test and paired t-test/Wilcoxon on reduction data."""
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
    
    return test_name, p_val, diff.mean()

def per_database_increase(df):
    """Calculates percentage increase for specific database columns."""
    print("\n" + "="*60)
    print("FUNCTIONAL GAINS: PERCENTAGE INCREASE PER DB")
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
    """Calculates the relative reduction percentage."""
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
# 3. PLOTTING FUNCTION
# ---------------------------------------------------------

def generate_plots(df, test_name, p_val, pct_increases, avg_total_improvement, OUTPUT_DIR):
    """Generates Figures 1 through 6."""
    metrics = ['GO_entries', 'COG_entries', 'KEGG_entries', 'EC_entries']

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
    
    # Only plot if we successfully added the Contigs column
    if 'Contigs' in df.columns and df['Contigs'].sum() > 0:
        df_sorted = df.sort_values(by='Contigs')
        
        plt.figure(figsize=(10, 6))
        
        # Plot Initial (Grey)
        plt.plot(df_sorted['Contigs'], df_sorted['Initial_Undesc_PCT'], 
                 marker='o', label='Initial', linestyle='-', color='grey', alpha=0.7)
        
        # Plot Enhanced (Red)
        plt.plot(df_sorted['Contigs'], df_sorted['Enhanced_Undesc_PCT'], 
                 marker='o', label='Enhanced', linestyle='-', color='red', alpha=0.7)
        
        plt.title("Impact of Assembly Quality on Annotations")
        plt.xlabel("Number of Contigs (calculated from GFFs)")
        plt.ylabel("Undescribed Coding Space (%)")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(f"{OUTPUT_DIR}/Figure3_Contigs_vs_Hypotheticals.png", dpi=300)
    else:
        print("Skipping Figure 3: No GFF directory provided or no contigs found.")

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
    # Recalculate if missing from df
    if 'Rel_Reduction_PCT' not in df.columns:
         raw_improvement = ((df['Initial_Undesc_PCT'] - df['Enhanced_Undesc_PCT']) / df['Initial_Undesc_PCT'].replace(0, 1)) * 100
         df['Rel_Reduction_PCT'] = raw_improvement.abs()

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

    # Figure 6: Percentage Increase per DB Bar Chart
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
        plt.ylim(0, db_df['Percent_Increase'].max() * 1.15) 
        plt.tight_layout()
        plt.savefig(f"{OUTPUT_DIR}/Figure6_DB_Percentage_Increase.png", dpi=300)
        print(f"-> Saved {OUTPUT_DIR}/Figure6_DB_Percentage_Increase.png")

    print(f"\nDone! Check the '{OUTPUT_DIR}' folder.")

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
        df = enrich_with_contigs(df, gff_dir_arg)
    else:
        print("3. [WARNING] No GFF directory detected (-g). Skipping contig counting.")

    # 2. Run Analysis
    test_name, p_val, _ = hypothetical_reduction(df)
    pct_increases = per_database_increase(df)
    avg_total_improvement = total_improvement(df)

    # 3. Create Output Directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 4. Generate Plots
    generate_plots(df, test_name, p_val, pct_increases, avg_total_improvement, output_dir)