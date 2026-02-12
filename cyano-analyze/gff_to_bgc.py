import os
import pandas as pd
import glob

def parse_gff_attributes(attribute_string):
    """Parses the GFF attribute column into a dictionary."""
    attrs = {}
    for item in attribute_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs

def parse_gff_file(filepath):
    """
    Parses a GFF file to extract BGC regions for antiSMASH and GECCO.
    Returns two lists of dictionaries: as_clusters, gecco_clusters
    """
    as_clusters = []
    gecco_clusters = []

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue

            tool_source = parts[1]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attributes = parse_gff_attributes(parts[8])

            # Filter for antiSMASH clusters
            if tool_source == 'antiSMASH' and feature_type == 'biosynthetic-gene-cluster':
                as_clusters.append({
                    'start': start,
                    'end': end,
                    'length': end - start,
                    'product': attributes.get('product', 'unknown'),
                    'id': attributes.get('ID', 'unknown')
                })
            
            # Filter for GECCO clusters
            elif tool_source == 'GECCO' and feature_type == 'biosynthetic-gene-cluster':
                gecco_clusters.append({
                    'start': start,
                    'end': end,
                    'length': end - start,
                    'product': attributes.get('product', 'unknown'),
                    'id': attributes.get('ID', 'unknown')
                })

    return as_clusters, gecco_clusters

def calculate_overlap_metrics(as_clusters, gecco_clusters):
    """Calculates overlap statistics between the two sets of clusters."""
    
    overlapping_as_count = 0
    total_overlap_bp = 0
    total_as_bp_covered = 0
    total_gecco_bp_covered = 0

    # Calculate overlaps (One AS can overlap multiple GECCOs and vice versa)
    for as_c in as_clusters:
        is_overlapping = False
        bp_covered_in_this_cluster = 0
        
        for gc_c in gecco_clusters:
            # Check for overlap: max(start1, start2) < min(end1, end2)
            overlap_start = max(as_c['start'], gc_c['start'])
            overlap_end = min(as_c['end'], gc_c['end'])
            
            if overlap_start < overlap_end:
                is_overlapping = True
                overlap_len = overlap_end - overlap_start
                total_overlap_bp += overlap_len
                
                # Simple accumulation of covered bases (approximation for multiple overlaps)
                bp_covered_in_this_cluster += overlap_len

        if is_overlapping:
            overlapping_as_count += 1
            # Cap the covered bp at the cluster length (in case of multi-overlaps)
            total_as_bp_covered += min(bp_covered_in_this_cluster, as_c['length'])

    # Reverse check to see how much of GECCO is covered by AntiSMASH
    for gc_c in gecco_clusters:
        bp_covered_in_this_cluster = 0
        for as_c in as_clusters:
            overlap_start = max(as_c['start'], gc_c['start'])
            overlap_end = min(as_c['end'], gc_c['end'])
            
            if overlap_start < overlap_end:
                overlap_len = overlap_end - overlap_start
                bp_covered_in_this_cluster += overlap_len
        
        total_gecco_bp_covered += min(bp_covered_in_this_cluster, gc_c['length'])

    return overlapping_as_count, total_overlap_bp, total_as_bp_covered, total_gecco_bp_covered

def analyze_directory(input_dir, output_file):
    results = []
    
    # Get all .gff or .gff3 files
    files = glob.glob(os.path.join(input_dir, "*.gff")) + glob.glob(os.path.join(input_dir, "*.gff3"))
    
    print(f"Found {len(files)} GFF files. Processing...")

    for filepath in files:
        filename = os.path.basename(filepath)
        locus_tag = os.path.splitext(filename)[0] # Assuming filename is the Locus Tag
        
        as_clusters, gecco_clusters = parse_gff_file(filepath)
        
        n_as = len(as_clusters)
        n_gecco = len(gecco_clusters)
        
        # Calculate metrics
        n_overlap, overlap_bp, as_bp_cov, gecco_bp_cov = calculate_overlap_metrics(as_clusters, gecco_clusters)
        
        # Derived Statistics
        total_as_length = sum(c['length'] for c in as_clusters) if n_as > 0 else 0
        total_gecco_length = sum(c['length'] for c in gecco_clusters) if n_gecco > 0 else 0
        
        # Coverage Percentages
        pct_as_covered = (as_bp_cov / total_as_length * 100) if total_as_length > 0 else 0
        pct_gecco_covered = (gecco_bp_cov / total_gecco_length * 100) if total_gecco_length > 0 else 0
        
        # Inter-tool Consensus (Common Loci / Total AntiSMASH Loci)
        consensus = (n_overlap / n_as) if n_as > 0 else 0
        
        # Overlap Efficiency (Dice Coefficient-like: 2*Overlap / (LenA + LenB))
        efficiency = (2 * overlap_bp / (total_as_length + total_gecco_length)) if (total_as_length + total_gecco_length) > 0 else 0
        
        # Consolidation Ratio (Fragmentation)
        consolidation = (n_as / n_gecco) if n_gecco > 0 else 0
        
        # Predominant Product (GECCO) - simple string concatenation for overview
        gecco_products = ", ".join(list(set([c['product'] for c in gecco_clusters])))
        
        results.append({
            'Locus Tag': locus_tag,
            'AntiSMASH Count': n_as,
            'GECCO Count': n_gecco,
            'Overlapping Clusters': n_overlap,
            'GECCO Products': gecco_products,
            'Total Overlap (bp)': overlap_bp,
            '% AntiSMASH Covered': round(pct_as_covered, 2),
            '% GECCO Covered': round(pct_gecco_covered, 2),
            'Inter-tool Consensus': round(consensus, 2),
            'Avg Overlap Efficiency': round(efficiency, 2),
            'Consolidation Ratio': round(consolidation, 2)
        })

    # Create DataFrame and Export
    df = pd.DataFrame(results)
    
    # Calculate Averages Row
    means = df.mean(numeric_only=True)
    means['Locus Tag'] = 'AVERAGE'
    df = pd.concat([df, pd.DataFrame([means])], ignore_index=True)

    df.to_excel(output_file, index=False)
    print(f"Analysis complete! Results saved to {output_file}")

# --- CONFIGURATION ---
# Άλλαξε αυτό το path στον φάκελο που έχεις τα GFF αρχεία σου
INPUT_DIRECTORY = "/home/angelos/Desktop/Thesis/test/full-mode/prochlorococaceae_gffs" 
OUTPUT_FILE = "Prochlorococaceae_BGC_Analysis_Results.xlsx"

if __name__ == "__main__":
    # Create a dummy folder and file for testing purposes if you run this script directly without setup
    # In real usage, just set INPUT_DIRECTORY correctly.
    if not os.path.exists(INPUT_DIRECTORY):
        print(f"Error: Directory {INPUT_DIRECTORY} not found.")
    else:
        analyze_directory(INPUT_DIRECTORY, OUTPUT_FILE)