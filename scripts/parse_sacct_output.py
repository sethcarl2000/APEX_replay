#!/usr/bin/env python3
"""
Modified from 'count_replayed_events' 

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

# the first line in a triplet of 3-SACCT output lines
line_main = re.compile(r'^apex_replay_(\d{4})_(\d{4})\|(\d+)_(\d+)\|(\w+)\|\|(\d{2}:\d{2}:\d{2}|\d{2}:\d{2}\.\d+)\|(\d{2}:\d{2}:\d{2}|\d{2}:\d{2}\.\d+)')

test_regex = re.compile(r'^apex_replay_(\d{4})_(\d{4})\|(\d+)_(\d+)\|(\w+)\|\|\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\|(\d{2}):(\d{2}):(\d{2})')


# the first line in a triplet of 3-SACCT output lines
line_batch = re.compile(r'^batch\|(\d{7})_(\d+).batch\|\w+\|(\d+K|\d+\.\d+M|\d+M)')


# output CSV format: 
# run-number-0 | run-number-1 | job-id | task-id | elapsed | total-cpu | max-RSS (MB) | job-state
FIELDS = ['run-number-0', 'run-number-1', 'job-id', 'task-id', 'elapsed', 'total-cpu', 'max-RSS-MB', 'job-state']

def parse_time(time_str):
    """Parse a time string in the format HH:MM:SS or MM:SS and return total seconds."""
    parts = time_str.split(':')
    if len(parts) == 3:
        hours, minutes, seconds = map(float, parts)
        return hours * 3600 + minutes * 60 + seconds
    elif len(parts) == 2:
        minutes, seconds = map(float, parts)
        return minutes * 60 + seconds
    else:
        raise ValueError(f"Invalid time format: {time_str}")

# ---------------------------------------------------------------
def parse_file(path): 
    """Yield one dict per block found in the file."""
    block = None

    with open(path, 'r', errors='replace') as infile:

        for line in infile:

            #print(f"line: {line.strip()}")  # debug print

            #m = test_regex.match(line)
            #if m:
            #    print(f"Matched test regex: {m.groups()}")  # debug print

            m = line_main.match(line)
            if m: 
                #print(f'Matched main line regex: {m.groups()}')  # debug print
                # A new job was found, start a new CSV line (flush existing line if we've already built one)
                if block is not None: 
                    #print(f"Yielding block.")  # debug print
                    yield block 
                
                # otherwise, fill out some info in this block 
                block = {f: None for f in FIELDS}

                block['run-number-0'] = int(m.group(1))
                block['run-number-1'] = int(m.group(2))
                block['job-id'] = int(m.group(3))
                block['task-id'] = int(m.group(4))
                block['job-state'] = m.group(5)
                elapsed = parse_time(m.group(6))
                total_cpu = parse_time(m.group(7))  # convert to seconds
                block['elapsed'] = elapsed
                block['total-cpu'] = total_cpu
                continue

            if block is None: 
                continue  # we haven't started a block yet

            # now, look for other lines
            m = line_batch.match(line)
            if m:
                #print(f'Matched batch line regex: {m.groups()}')  # debug print
                
                # first, check to make sure the job & task id are the same
                new_job_id = int(m.group(1))
                new_task_id = int(m.group(2))
                if new_job_id != block['job-id'] or new_task_id != block['task-id']:
                    print(f"Warning: job/task mismatch in {path}: "
                          f"expected {block['job-id']}/{block['task-id']}, "
                          f"found {new_job_id}/{new_task_id}")
                    continue
                
                max_RSS = m.group(3)
                if max_RSS[-1] == 'K':
                    block['max-RSS-MB'] = float(max_RSS[:-1]) / 1024.0  # convert KB to MB
                elif max_RSS[-1] == 'M':
                    block['max-RSS-MB'] = float(max_RSS[:-1])  # already in MB
                else:
                    print(f"Warning: unrecognized max-RSS format in {path}: {max_RSS}")
                    block['max-RSS-MB'] = None

                continue 
        
        if block is not None:
            yield block

# ---------------------------------------------------------------       
0

def main():
    ap = argparse.ArgumentParser(description=__doc__)

    ap.add_argument('input_file', type=Path, help='path to the input file')

    ap.add_argument('-o', '--output', type=Path, default=Path('misc/jobs.csv'),
                    help='output CSV path (default: misc/jobs.csv)')
    args = ap.parse_args()

    n_jobs = 0
    n_discarded = 0 
    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=FIELDS)
        writer.writeheader()

        for block in parse_file(args.input_file):
            if any(block[f] is None for f in FIELDS):
                print(f"Warning: discarded incomplete job block: {block}", file=sys.stderr)
                n_discarded += 1
                continue
            writer.writerow(block)
            n_jobs += 1

    print(f"Found {n_jobs} job(s); wrote {args.output}")

    if n_discarded:
        print(f"Warning: discarded {n_discarded} incomplete job(s) "
              f"(missing the 'batch' job info from slurm for this job)", file=sys.stderr)




if __name__ == '__main__':
    main()
