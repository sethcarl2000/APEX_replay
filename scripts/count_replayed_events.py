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
RE_DECODE = re.compile(r'\[([0-9:.]+)\]\s*Info\s*from\s*<execute_specific_run>:\s*Processing\s*decode\.')
RE_DONE_DECODE = re.compile(r'^.*\[([0-9:.]+)\]\s*Info\s*from\s*<execute_specific_run>:\s*Done\s*with\s*decode')

RUN_RE = re.compile(r'--output\s*file\s*/volatile/halla/apex/full_replay/production/replay\.(\d+)\.rawfile-(\d+)\.segment-(\d+)\.root')

RE_REPLAY_START = re.compile(r'^.*\[([0-9:.]+)\]\s*Info\s*from\s*<Process>:\s*In\s*funciton\s*call\s*body\.')
RE_REPLAY_END = re.compile(r'^.*\[([0-9:.]+)\]\s*Info\s*from\s*<Process>:\s*(?:Acceptable\s*replay\s*completed)')

EVENTS_RE = re.compile(r'^.*Processing\s+(\d+)\s+events') # Processing\s+(\d+)\s+events')

# Ordered list of (csv_column_name, regex). The ORDER here must match 
# the order the lines appear within a block. Duplicate patterns are
# fine: the first occurrence in the block fills the first entry with
# that pattern, the second occurrence fills the next, etc.
VAR_SPECS = [
    ('one_s2_both_arms', re.compile(r'\s+1 S2 hit in each arm:\s+(\d+)')),
    ('s2_coinc', re.compile(r'\s+1 coinc with LHRS & RHRS:\s+(\d+)')),
    ('one_group_RHRS', re.compile(r'\s+1 hit group\s+(\d+)')),
    ('one_pair_RHRS',  re.compile(r'\s+1 hit group pair\s+(\d+)')),
    ('one_raw_track_RHRS', re.compile(r'\s+1 raw track\s+(\d+)')),
    ('one_track_RHRS', re.compile(r'\s+1 refined track\s+(\d+)')),
    ('one_group_LHRS', re.compile(r'\s+1 hit group\s+(\d+)')),
    ('one_pair_LHRS',  re.compile(r'\s+1 hit group pair\s+(\d+)')),
    ('one_raw_track_LHRS', re.compile(r'\s+1 raw track\s+(\d+)')),
    ('one_track_LHRS', re.compile(r'\s+1 refined track\s+(\d+)')),  # same regex as var_pattern_2
]

# ---------------------------------------------------------------

VAR_NAMES = [name for name, _ in VAR_SPECS]
FIELDS = ['run', 'rawfile_num', 'segment_num', 'events', 'decode_time', 'replay_time'] + VAR_NAMES


def parse_timestamp(tstamp: str):

    timestamp_re = re.compile(r'(\d+):(\d+):([0-9.]+)')

    if (m := timestamp_re.match(tstamp)) is not None: 
        seconds = float(m.group(1).strip())*3600.
        seconds = float(m.group(2).strip())*60. + seconds
        seconds = float(m.group(3).strip()) + seconds
        return seconds
    else: 
        return None 


def parse_file(path):

    """Yield one dict per block found in the file."""
    block = None

    start_decode = None
    start_replay = None 

    with open(path, 'r', errors='replace') as fh:

        print(f"processing file: {path}") 

        for line in fh:

            m = RE_DECODE.match(line)
            if m: 
                # A new run line starts a new block; flush any previous one
                if block is not None:
                    yield block
                block = {f: None for f in FIELDS}   

                if (start_decode := parse_timestamp(m.group(1).strip())) is None: 
                    raise ValueError(f"unable to parse timestamp: {m.group(1).strip()}")

                start_replay = None 
                block['decode_time'] = "null"
                continue
            
            if block is None:
                continue  # not inside a block yet

            m = RE_DONE_DECODE.match(line)
            if m:

                if (end_decode := parse_timestamp(m.group(1).strip())) is None: 
                    ValueError(f"unable to parse timestamp: {m.group(1).strip()}")

                block['decode_time'] = f"{end_decode - start_decode:.4f}"   
                continue 

            if (m := RE_REPLAY_START.match(line)) is not None: 

                if (start_replay := parse_timestamp(m.group(1).strip())) is None: 
                    ValueError(f"unable to parse timestamp: {m.group(1).strip()}")

                continue 

            if (m := RE_REPLAY_END.match(line)) is not None and start_replay is not None: 
            
                if (end_replay := parse_timestamp(m.group(1).strip())) is None: 
                    ValueError(f"unable to parse timestamp: {m.group(1).strip()}")

                block['replay_time'] = f"{end_replay - start_replay:.4f}"
                continue 
            

            m = RUN_RE.match(line)
            if m:
                block['run'] = int(m.group(1))
                block['rawfile_num'] = int(m.group(2))
                block['segment_num'] = int(m.group(3))
                continue


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
