import argparse
from . import bgc
from . import plots
from . import prokka_update

def main():
    parser = argparse.ArgumentParser(description="Genome Analysis Package")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ---------------------------------------------------------
    # COMMAND: stats (CSV Analysis)
    # ---------------------------------------------------------
    p_stats = subparsers.add_parser("stats", help="Generate stats/plots from CSV data")
    p_stats.add_argument("-i", "--input", required=True, help="Input CSV file")
    p_stats.add_argument("-o", "--output", default="plots", help="Output directory")
    p_stats.add_argument("-g", "--gff-dir", help="Directory containing GFF files")

    # ---------------------------------------------------------
    # COMMAND: bgc (GFF/Contig Analysis)
    # ---------------------------------------------------------
    p_bgc = subparsers.add_parser("bgc", help="Analyze BGCs and Contigs from a directory of GFF files")
    
    # CHANGED: Input is now strictly a directory
    p_bgc.add_argument("-d", "--dir", required=True, help="Directory containing .gff or .gff3 files")
    p_bgc.add_argument("-o", "--output", default="bgc_results.xlsx", help="Output Excel file")
    p_bgc.add_argument("--skip-plots", action="store_true", help="Skip statistical plots")

    # ---------------------------------------------------------
    # COMMAND: prokka (Update CSV with Prokka Metrics)
    # ---------------------------------------------------------
    p_prokka = subparsers.add_parser("prokka", help="Add Prokka metrics to an existing summary CSV")
    p_prokka.add_argument("--csv_file", required=True, help="Path to the existing summary CSV")
    p_prokka.add_argument("--prokka_gff_dir", required=True, help="Directory containing Prokka GFF files")
    p_prokka.add_argument("--prokka_gbk_dir", required=True, help="Directory containing Prokka GBK files")

    # ---------------------------------------------------------
    # PARSING & EXECUTION
    # ---------------------------------------------------------
    args = parser.parse_args()

    if args.command == "stats":
        plots.cmd_plot(args)

    elif args.command == "bgc":
        bgc.cmd_plot(args)

    elif args.command == "prokka":
        # Call the main function from your new module
        prokka_update.update_all_genomes(
            args.csv_file,
            args.prokka_gff_dir,
            args.prokka_gbk_dir
        )


if __name__ == "__main__":
    main()