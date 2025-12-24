#!/usr/bin/env python3
"""Cross-platform cleanup script for opm_align results."""

import os
import sys
import shutil
import glob


RESULTS_DIR = "results"


def print_usage():
    """Print usage information."""
    print(f"Usage: {sys.argv[0]} [fasta|log|db|all]")
    print("Remove fasta, log or all files from results directory")
    print("  fasta: remove fasta files")
    print("  db:    remove db files")
    print("  log:   remove log files")
    print("  all:   remove all files")


def remove_file(filepath):
    """Remove a file if it exists."""
    if os.path.isfile(filepath):
        os.remove(filepath)
        print(f"  Removed: {filepath}")
    else:
        print(f"  Not found: {filepath}")


def remove_files_pattern(pattern):
    """Remove files matching a glob pattern."""
    files = glob.glob(pattern)
    for f in files:
        if os.path.isfile(f):
            os.remove(f)
            print(f"  Removed: {f}")


def cleanup_fasta():
    """Remove fasta files."""
    print("Removing fasta files")
    remove_file(os.path.join(RESULTS_DIR, "db_fasta.txt"))
    remove_file(os.path.join(RESULTS_DIR, "new_fasta.txt"))


def cleanup_log():
    """Remove log files."""
    print("Removing log files")
    remove_file(os.path.join(RESULTS_DIR, "find_best_match.log"))
    remove_file(os.path.join(RESULTS_DIR, "warnings.log"))


def cleanup_db():
    """Remove all db files."""
    print("Removing all db files")
    blast_db_path = os.path.join(RESULTS_DIR, "blast_db")
    if os.path.exists(blast_db_path) or glob.glob(os.path.join(RESULTS_DIR, "db*")):
        remove_files_pattern(os.path.join(RESULTS_DIR, "db*"))
    else:
        print("  blast_db directory does not exist")


def cleanup_all():
    """Remove all files in results directory."""
    print("Removing all files")
    if os.path.isdir(RESULTS_DIR):
        for item in os.listdir(RESULTS_DIR):
            item_path = os.path.join(RESULTS_DIR, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
                print(f"  Removed file: {item_path}")
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
                print(f"  Removed directory: {item_path}")
    else:
        print("  results directory does not exist")


def main():
    """Main entry point."""
    if len(sys.argv) != 2:
        print_usage()
        sys.exit(1)

    command = sys.argv[1].lower()

    if command == "fasta":
        cleanup_fasta()
    elif command == "log":
        cleanup_log()
    elif command == "db":
        cleanup_db()
    elif command == "all":
        cleanup_all()
    else:
        print_usage()
        sys.exit(1)


if __name__ == "__main__":
    main()
