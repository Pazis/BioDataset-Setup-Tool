import argparse
from . import bgc
from . import plots
from . import prokka_update
import sys
import os
import logging

class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass

def setup_tool_logging(tool_name, output_dir):
    """Sets up a log file in the output directory and redirects stdout/stderr."""
    # Ensure output directory exists so we can write the log
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    log_file = os.path.join(output_dir, f"{tool_name}.log")

    # Configure the Root Logger
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'), # Write to file
            logging.StreamHandler(sys.stdout)        # Keep writing to console
        ],
        force=True
    )
    
    # Redirect standard print() calls to the logger
    sys.stdout = StreamToLogger(logging.getLogger(tool_name), logging.INFO)
    sys.stderr = StreamToLogger(logging.getLogger(f"{tool_name}_ERROR"), logging.ERROR)
    
    logging.info(f"--- Starting {tool_name} Analysis ---")
    logging.info(f"Log file location: {log_file}")


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
        setup_tool_logging("stats", args.output)
        plots.cmd_plot(args)

    elif args.command == "bgc":
        out_dir = args.output if os.path.isdir(args.output) else os.path.dirname(args.output)
        if not out_dir: out_dir = "." # Default to current dir if no path provided
        
        setup_tool_logging("bgc", out_dir)
        bgc.cmd_plot(args)

    elif args.command == "prokka":
        # Call the main function from your new module
        out_dir = os.path.dirname(args.csv_file)
        if not out_dir: out_dir = "."
        
        setup_tool_logging("prokka", out_dir)
        prokka_update.update_all_genomes(
            args.csv_file,
            args.prokka_gff_dir,
            args.prokka_gbk_dir
        )

if __name__ == "__main__":
    main()