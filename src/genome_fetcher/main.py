import argparse
import time
import shutil
import zipfile
import pandas as pd
import requests
import xml.etree.ElementTree as ET
from pathlib import Path
from tqdm import tqdm

# Ρυθμίσεις
TAXON = "cyanobacteria"  
OUTPUT_RAW = "cyanobacteria_raw.csv"
OUTPUT_FILTERED = "cyanobacteria_selected.xlsx"
ZIP_FILE = "genomes_dataset.zip"
FINAL_DATASET_DIR = "Cyano_Benchmark_Dataset" # Εδώ θα μπουν τα τακτοποιημένα

# API Endpoints
API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ENTREZ_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def get_lineage_batch(taxids):
    """Ανάκτηση taxonomy από NCBI Entrez"""
    results = {}
    batch_size = 100
    print(f"🧬 Ανάκτηση ταξινόμησης για {len(taxids)} μοναδικά TaxIDs...")
    
    for i in tqdm(range(0, len(taxids), batch_size), desc="Taxonomy Batches"):
        batch = taxids[i:i + batch_size]
        id_str = ",".join(map(str, batch))
        params = {"db": "taxonomy", "id": id_str, "retmode": "xml"}
        try:
            r = requests.post(ENTREZ_URL, data=params)
            r.raise_for_status()
            root = ET.fromstring(r.text)
            for taxon in root.findall("Taxon"):
                tid = taxon.find("TaxId").text
                family, order = "N/A", "N/A"
                lineage_ex = taxon.find("LineageEx")
                if lineage_ex is not None:
                    for entry in lineage_ex.findall("Taxon"):
                        rank = entry.find("Rank").text
                        name = entry.find("ScientificName").text
                        if rank == "family": family = name
                        elif rank == "order": order = name
                results[tid] = {"family": family, "order": order}
            time.sleep(0.35)
        except Exception as e:
            print(f"⚠️ Error fetching batch: {e}")
    return results

def cmd_fetch(args):
    print(f"🚀 Ξεκινάω fetch metadata για: {TAXON}...")
    url = f"{API_BASE}/genome/taxon/1117/dataset_report"
    params = {"page_size": 1000}
    all_genomes = []
    next_token = None
    page_num = 1
    
    try:
        while True:
            if next_token: params["page_token"] = next_token
            print(f"⏳ Λήψη σελίδας {page_num}...")
            response = requests.get(url, params=params)
            response.raise_for_status()
            data_json = response.json()
            reports = data_json.get('reports', [])
            if not reports and page_num == 1: return

            for record in reports:
                asm_info = record.get('assembly_info', {})
                org = record.get('organism', {})
                acc = record.get('accession') or record.get('current_accession')
                if not acc: continue
                all_genomes.append({
                    "accession": acc,
                    "taxid": str(org.get('tax_id')),
                    "organism": org.get('organism_name'),
                    "assembly_level": asm_info.get('assembly_level', 'N/A'),
                    "submission_date": asm_info.get('release_date', '1900-01-01'),
                    "family": "N/A", "order": "N/A"
                })
            next_token = data_json.get('next_page_token')
            if not next_token: break
            page_num += 1
    except Exception as e:
        print(f"❌ Error: {e}")
        return

    unique_taxids = list(set(g['taxid'] for g in all_genomes if g['taxid']))
    taxonomy_map = get_lineage_batch(unique_taxids)
    for genome in all_genomes:
        tid = genome['taxid']
        if tid in taxonomy_map:
            genome['family'] = taxonomy_map[tid]['family']
            genome['order'] = taxonomy_map[tid]['order']

    df = pd.DataFrame(all_genomes)
    df.to_csv(OUTPUT_RAW, index=False)
    print(f"💾 Αποθηκεύτηκαν στο {OUTPUT_RAW}")

def cmd_filter(args):
    print(f"🔍 Φιλτράρισμα: Κρατάμε {args.number} ανά Οικογένεια...")
    if not Path(OUTPUT_RAW).exists(): return
    df = pd.read_csv(OUTPUT_RAW)
    df = df[df['family'] != 'N/A']
    
    quality_map = {'Complete Genome': 3, 'Chromosome': 2, 'Scaffold': 1, 'Contig': 0}
    df['quality_score'] = df['assembly_level'].map(quality_map).fillna(0)
    df['submission_date'] = pd.to_datetime(df['submission_date'])
    
    df_sorted = df.sort_values(by=['quality_score', 'submission_date'], ascending=[False, False])
    selected_df = df_sorted.groupby('family').head(args.number)
    
    print(selected_df['family'].value_counts().head(5))
    selected_df.to_excel(OUTPUT_FILTERED, index=False)
    
    with open("accession_list.txt", "w") as f:
        for acc in selected_df['accession']: f.write(f"{acc}\n")
    print(f"✅ Λίστα έτοιμη: {OUTPUT_FILTERED}")

def cmd_download(args):
    print("⬇️ Έναρξη Download...")
    if not Path("accession_list.txt").exists(): return
    with open("accession_list.txt", "r") as f:
        accessions = [line.strip() for line in f if line.strip()]
    
    url = f"{API_BASE}/genome/download"
    payload = {
        "accessions": accessions,
        "include_annotation_type": ["GENOME_FASTA", "GENOME_GFF", "PROT_FASTA"]
    }
    
    try:
        print(f"⏳ Downloading ZIP ({len(accessions)} genomes)...")
        with requests.post(url, json=payload, stream=True) as r:
            r.raise_for_status()
            with open(ZIP_FILE, 'wb') as f:
                total_size = int(r.headers.get('content-length', 0))
                with tqdm(total=total_size, unit='B', unit_scale=True, desc="Download") as pbar:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
        print(f"✅ Download ολοκληρώθηκε: {ZIP_FILE}")
    except Exception as e:
        print(f"❌ Error: {e}")

def cmd_organize(args):
    """
    Αυτή η εντολή κάνει unzip και τακτοποιεί τα αρχεία ΜΟΝΟ σε:
    - high_quality (Complete Genome, Chromosome)
    - low_quality (Scaffold, Contig)
    Χωρίς υποφακέλους οικογενειών.
    """
    print("📂 Οργάνωση Dataset σε High/Low Quality...")
    
    # 1. Έλεγχοι αρχείων
    if not Path(ZIP_FILE).exists():
        print(f"❌ Δεν βρέθηκε το {ZIP_FILE}. Τρέξε πρώτα 'download'.")
        return
    if not Path(OUTPUT_FILTERED).exists():
        print(f"❌ Δεν βρέθηκε το Excel {OUTPUT_FILTERED}. Τρέξε πρώτα 'filter'.")
        return

    # 2. Φόρτωση Mapping
    print("📖 Διάβασμα Excel για αντιστοίχιση...")
    df = pd.read_excel(OUTPUT_FILTERED)
    
    # Χρειαζόμαστε μόνο το Assembly Level τώρα
    # Το accession πρέπει να υπάρχει στο excel για να το συμπεριλάβουμε
    acc_to_level = dict(zip(df['accession'], df['assembly_level']))

    # 3. Unzip σε προσωρινό φάκελο
    temp_dir = Path("temp_ncbi_unzipped")
    if temp_dir.exists(): shutil.rmtree(temp_dir)
    
    print("⏳ Unzipping (μπορεί να πάρει λίγο)...")
    with zipfile.ZipFile(ZIP_FILE, 'r') as zip_ref:
        zip_ref.extractall(temp_dir)

    # 4. Δημιουργία Τελικής Δομής
    final_dir = Path(FINAL_DATASET_DIR)
    if final_dir.exists(): shutil.rmtree(final_dir)
    final_dir.mkdir()
    
    # Δημιουργία υποφακέλων ποιότητας μόνο
    hq_dir = final_dir / "high_quality"
    lq_dir = final_dir / "low_quality"
    
    hq_dir.mkdir()
    lq_dir.mkdir()

    print("🚀 Μετακίνηση και οργάνωση αρχείων...")
    
    data_path = temp_dir / "ncbi_dataset" / "data"
    count_moved = 0
    
    # Iteration σε όλους τους φακέλους
    for genome_dir in tqdm(list(data_path.iterdir())):
        if not genome_dir.is_dir(): continue
        
        accession = genome_dir.name
        
        # Αν δεν είναι στη λίστα μας (filtered), το αγνοούμε
        level = acc_to_level.get(accession)
        if not level:
            continue

        # Καθορισμός Quality Folder
        level_lower = str(level).lower()
        if 'complete' in level_lower or 'chromosome' in level_lower:
            target_folder = hq_dir
        else:
            target_folder = lq_dir

        # Αντιγραφή αρχείων απευθείας στον φάκελο ποιότητας
        for file_path in genome_dir.glob("*"):
            if file_path.is_file():
                shutil.copy(file_path, target_folder / file_path.name)
        
        count_moved += 1

    # 5. Καθαρισμός
    print("🧹 Καθαρισμός προσωρινών αρχείων...")
    shutil.rmtree(temp_dir)
    
    print(f"\n✅ Ολοκληρώθηκε! Οργανώθηκαν {count_moved} γονιδιώματα.")
    print(f"📂 Δομή φακέλων: {FINAL_DATASET_DIR} -> [high_quality / low_quality]")

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("fetch")
    p_filter = subparsers.add_parser("filter")
    p_filter.add_argument("-n", "--number", type=int, default=10) # Default 10 για benchmark
    subparsers.add_parser("download")
    subparsers.add_parser("organize") # Νέα εντολή

    args = parser.parse_args()

    if args.command == "fetch":
        cmd_fetch(args)
    elif args.command == "filter":
        cmd_filter(args)
    elif args.command == "download":
        cmd_download(args)
    elif args.command == "organize":
        cmd_organize(args)

if __name__ == "__main__":
    main()