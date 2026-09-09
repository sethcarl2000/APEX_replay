#!/usr/bin/env python3
"""
collect sizes of all rawfiles, and place them into a csv

Usage:
    python3 scan_blocks.py /path/to/directory [-o output.csv]
"""

import argparse
import csv
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------
# Configuration: edit these to match your actual variable names
# ---------------------------------------------------------------

# the first line in a triplet of 3-SACCT output lines
line_patt = re.compile(r'^(\d+)\s+/cache/halla/apex/raw/apex_(\d{4}).dat.(\d+)')

# output CSV format: 
# run-number | rawfile-number | size-MB
FIELDS = ['run-number', 'rawfile-number', 'size-MB']

# ---------------------------------------------------------------       
0

def main():
    ap = argparse.ArgumentParser(description=__doc__)

    ap.add_argument('input_file', type=Path, help='path to the input file')

    ap.add_argument('-o', '--output', type=Path, default=Path('misc/jobs.csv'),
                    help='output CSV path (default: misc/jobs.csv)')
    args = ap.parse_args()

    n_files = 0
    n_discarded = 0 
    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS)
        writer.writeheader()

        with open(args.input_file, 'r') as infile: 
            for line in infile: 
                m = line_patt.match(line)
                if m: 
                    run_number = int(m.group(2))
                    rawfile_number = int(m.group(3))
                    size_MB = float(m.group(1)) / (1024.0)  # convert bytes to MB
                    block = {
                        'run-number': run_number,
                        'rawfile-number': rawfile_number,
                        'size-MB': size_MB
                    }
                    writer.writerow(block)
                    n_files += 1
                else: 
                    print(f"Warning: line did not match expected format: {line.strip()}", file=sys.stderr)
                    n_discarded += 1

    print(f"Found {n_files} file(s); wrote {args.output}")

    if n_discarded:
        print(f"Warning: discarded {n_discarded} incomplete lines", file=sys.stderr)




if __name__ == '__main__':
    main()
