#!/usr/bin/env python3
"""
Written by Claude Fable 5 7/9/2026

Scan text files in a directory for analysis blocks, extract figures,
write one CSV row per block, and print cumulative sums.

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

# Matches "--run nubmer 1234" (tolerates both 'nubmer' and 'number')
RUN_RE = re.compile(r'^^\s*Processing run/rawfile/segment:\s+(\d+)')

EVENTS_RE = re.compile(r'^\s*Processing\s+(\d+)\s+events')

# Ordered list of (csv_column_name, regex). The ORDER here must match
# the order the lines appear within a block. Duplicate patterns are
# fine: the first occurrence in the block fills the first entry with
# that pattern, the second occurrence fills the next, etc.
VAR_SPECS = [
    ('one_s2_both_arms', re.compile(r'^\s*1 S2 hit in each arm:\s+(\d+)')),
    ('s2_coinc', re.compile(r'^\s*1 coinc with LHRS & RHRS:\s+(\d+)')),
    ('one_track_RHRS', re.compile(r'^\s*1 refined track\s+(\d+)')),
    ('one_track_LHRS', re.compile(r'^\s*1 refined track\s+(\d+)')),  # same regex as var_pattern_2
]

# ---------------------------------------------------------------

VAR_NAMES = [name for name, _ in VAR_SPECS]
FIELDS = ['file', 'run', 'events'] + VAR_NAMES


def parse_file(path):
    """Yield one dict per block found in the file."""
    block = None
    with open(path, 'r', errors='replace') as fh:

        print("opened path: ", path) 

        for line in fh:
            m = RUN_RE.match(line)
            if m:
                # A new run line starts a new block; flush any previous one
                if block is not None:
                    yield block
                block = {f: None for f in FIELDS}
                block['file'] = path.name
                block['run'] = int(m.group(1))
                continue

            if block is None:
                continue  # not inside a block yet

            m = EVENTS_RE.match(line)
            if m:
                block['events'] = int(m.group(1))
                continue

            # Assign to the FIRST still-empty column whose pattern matches.
            # This is what disambiguates identical patterns by position.
            for name, pat in VAR_SPECS:
                if block[name] is not None:
                    continue
                m = pat.match(line)
                if m:
                    block[name] = int(m.group(1))
                    break

    if block is not None:
        yield block


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('directory', type=Path, help='directory containing text files')
    ap.add_argument('-o', '--output', type=Path, default=Path('blocks.csv'),
                    help='output CSV path (default: blocks.csv)')
    ap.add_argument('-g', '--glob', default='*',
                    help="filename glob to match, e.g. '*.log' (default: all files)")
    args = ap.parse_args()

    if not args.directory.is_dir():
        sys.exit(f"Error: {args.directory} is not a directory")

    sum_fields = ['events'] + VAR_NAMES
    totals = {f: 0 for f in sum_fields}
    n_blocks = 0
    n_discarded = 0

    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS)
        writer.writeheader()

        for path in sorted(args.directory.glob(args.glob)):
            if not path.is_file():
                continue
            for block in parse_file(path):
                # --- NEW: discard blocks with any missing value ---
                if any(block[f] is None for f in FIELDS):
                    n_discarded += 1
                    continue
                # --------------------------------------------------
                writer.writerow(block)
                n_blocks += 1
                for f in sum_fields:
                    if block[f] is not None:
                        totals[f] += block[f]
                    else:
                        missing += 1

    print(f"Found {n_blocks} block(s); wrote {args.output}")

    if n_discarded:
        print(f"Warning: discarded {n_discarded} incomplete block(s) "
              f"(missing run, event count, or one of the variables)", file=sys.stderr)



    print("\nCumulative sums:")
    width = max(len(f) for f in sum_fields)
    for f in sum_fields:
        print(f"  {f:<{width}}  {totals[f]:>15,}")


if __name__ == '__main__':
    main()
