#!/usr/bin/env python3
"""
slurm_hist.py — Graphical histogram wrapper for SLURM's sacct command.

Usage examples
--------------
# Histogram of Elapsed and MaxRSS for all jobs by user 'alice' since Jan 1:
python3 slurm_hist.py --user alice --starttime 2025-01-01 --fields Elapsed MaxRSS

# Filter by job-name wildcard, e.g. 'sim_run_*':
python3 slurm_hist.py --user alice --starttime 2025-01-01 \
    --jobname "sim_run_*" --fields Elapsed MaxRSS CPUTime

# Save figures instead of displaying them:
python3 slurm_hist.py --user alice --starttime 2025-01-01 \
    --fields Elapsed MaxRSS --outdir ./plots

# Change number of histogram bins:
python3 slurm_hist.py --user alice --starttime 2025-01-01 \
    --fields Elapsed --bins 30
"""

import argparse
import subprocess
import sys
import re
import os
from pathlib import Path
from fnmatch import fnmatch

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── sacct field helpers ────────────────────────────────────────────────────────

# Fields whose values are HH:MM:SS (or D-HH:MM:SS) — convert to seconds.
TIME_FIELDS = {
    "Elapsed", "CPUTime", "CPUTimeRAW", "AveCPU", "MinCPU", "SystemCPU",
    "UserCPU", "Timelimit", "Reserved", "ReqCPUFreqMax", "ReqCPUFreqMin",
}

# Fields whose values carry a memory suffix (K / M / G / T) — convert to MiB.
MEM_FIELDS = {
    "MaxRSS", "MaxVMSize", "AveRSS", "AveVMSize", "MaxDiskRead",
    "MaxDiskWrite", "AveDiskRead", "AveDiskWrite", "ReqMem",
}

# Human-readable unit labels after conversion
UNIT_LABEL = {
    **{f: "seconds" for f in TIME_FIELDS},
    **{f: "MiB"    for f in MEM_FIELDS},
}


def parse_time(value: str) -> float:
    """Convert [D-]HH:MM:SS to total seconds. Returns None if unparseable."""
    value = value.strip()
    if not value or value in ("", "INVALID", "Unknown", "None"):
        return None
    days = 0
    if "-" in value:
        d, value = value.split("-", 1)
        days = int(d)
    parts = value.split(":")
    try:
        if len(parts) == 3:
            h, m, s = parts
            return days * 86400 + int(h) * 3600 + int(m) * 60 + float(s)
        if len(parts) == 2:
            m, s = parts
            return days * 86400 + int(m) * 60 + float(s)
    except ValueError:
        pass
    return None


def parse_memory(value: str) -> float:
    """Convert sacct memory string (e.g. '4096K', '2.5G') to MiB. Returns None if unparseable."""
    
    value = value.strip()
    if not value or value in ("", "0", "INVALID", "Unknown", "None"):
        return None
    multipliers = {"K": 1 / 1024, "M": 1.0, "G": 1024.0, "T": 1024.0 ** 2}
    suffix = value[-1]# .upper()
    
    if suffix in multipliers:
        try:
            return float(value[:-1]) * multipliers[suffix]
        except ValueError:
            pass
    try:
        # Assume bytes if no suffix
        return float(value) / (1024 ** 2)
    except ValueError:
        return None


def parse_numeric(value: str) -> float:
    """Generic numeric parse; returns None for empty/sentinel values."""
    value = value.strip()
    if not value or value in ("", "INVALID", "Unknown", "None"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def convert_field(field: str, raw: str) -> float:
    """Dispatch a raw sacct string to the correct converter for the given field."""
    if field in TIME_FIELDS:
        #print(f"field {field} is Time")
        return parse_time(raw)
    if field in MEM_FIELDS:
        #print(f"field {field} is Mem")
        return parse_memory(raw)
    return parse_numeric(raw)


# ── sacct invocation ───────────────────────────────────────────────────────────

def build_sacct_cmd(args: argparse.Namespace, fields: list[str]) -> list[str]:
    """Construct the sacct command list from parsed CLI arguments."""
    # Always include JobName and State so we can filter on them
    format_fields = ["JobName", "State"] + fields
    format_str = ",".join(format_fields)

    cmd = [
        "sacct",
        "--parsable2",          # pipe-delimited, no trailing |
        "--noheader",           # strip off the header
        f"--format={format_str}",
        "--state=cd",           # completed jobs only (override with --state)
        "--endtime=now",        # by defualt, show all jobs up to 'now'. 
    ]

    if args.user:
        cmd += [f"--user={args.user}"]
    if args.starttime:
        cmd += [f"--starttime={args.starttime}"]
    if args.endtime:
        # Remove our default --state=CD if the user supplies their own
        cmd = [c for c in cmd if not c.startswith("--endtime=")]
        cmd += [f"--endtime={args.endtime}"]
    if args.account:
        cmd += [f"--account={args.account}"]
    if args.partition:
        cmd += [f"--partition={args.partition}"]
    if args.state:
        # Remove our default --state=CD if the user supplies their own
        cmd = [c for c in cmd if not c.startswith("--state=")]
        cmd += [f"--state={args.state}"]
    if args.jobs:
        cmd += [f"--jobs={','.join(args.jobs)}"]

    return cmd


def run_sacct(cmd: list[str]) -> list[list[str]]:
    """Run sacct and return rows as lists of strings."""
    try:
        print(f"sacct command: {cmd}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        sys.exit(
            "ERROR: 'sacct' not found. Are you running this on a SLURM login/compute node?"
        )
    except subprocess.CalledProcessError as exc:
        sys.exit(f"ERROR: sacct exited with code {exc.returncode}:\n{exc.stderr}")

    rows = []
    for line in result.stdout.splitlines():
        parts = line.split("|")
        if parts:
            rows.append(parts)
    return rows


# ── data parsing ───────────────────────────────────────────────────────────────

def collect_data(
    rows: list[list[str]],
    fields: list[str],
    jobname_pattern: str,
) -> dict[str, list[float]]:
    """
    Parse sacct rows into per-field lists of floats.

    sacct returns multiple rows per job (one per job-step, e.g. 'batch', 'extern').
    We keep only the base job row (JobName matches directly, not 'batch'/'extern').
    Optionally filter by a wildcard job-name pattern.
    """
    # Column indices: [JobName, State, field0, field1, ...]
    data: dict[str, list[float]] = {f: [] for f in fields}
    field_indices = {f: i + 2 for i, f in enumerate(fields)}

    # each job is a 'triplet' of 3 jobs
    n_results = int(len(rows)/3)

    if len(rows) % 3 != 0:
        sys.exit(f"ERROR: number of rows ({len(rows)}) is not divisible by 3") 

    for i_result in range(0, n_results):
        
        job_name_row = rows[(3*i_result) + 0] # the first row, which gives us the job name
        data_row     = rows[(3*i_result) + 1] 

        job_name = job_name_row[0]
        state    = job_name_row[1]
        
        # Wildcard filter
        if jobname_pattern and not fnmatch(job_name, jobname_pattern):
            continue

        if len(data_row) < 2 + len(fields):
            continue    
        
        for field in fields:
            idx = field_indices[field]
            if idx < len(data_row):
                val = convert_field(field, data_row[idx])
                if val is not None:
                    data[field].append(val)

    return data


# ── plotting ───────────────────────────────────────────────────────────────────

# A clean, terminal/data-science palette
PALETTE = [
    "#4C9BE8", "#E8834C", "#5DBF76", "#C96DD8",
    "#E8C14C", "#4CCFE8", "#E84C6D", "#A8E84C",
]

def format_axis_label(field: str) -> str:
    unit = UNIT_LABEL.get(field, "")
    return f"{field} ({unit})" if unit else field


def plot_histogram(
    field: str,
    values: list[float],
    bins: int,
    ax: plt.Axes,
    color: str,
    jobname_pattern: str, # str | None,
) -> None:
    arr = np.array(values)

    n, bin_edges, patches = ax.hist(
        arr,
        bins=bins,
        color=color,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.88,
    )

    # Overlay a KDE for smooth shape hint (only when enough data)
    if 0: # len(arr) >= 10:
        from scipy.stats import gaussian_kde  # optional; graceful fallback below
        kde_x = np.linspace(arr.min(), arr.max(), 300)
        kde_y = gaussian_kde(arr)(kde_x)
        # Scale KDE to histogram height
        bin_width = bin_edges[1] - bin_edges[0]
        kde_y_scaled = kde_y * len(arr) * bin_width
        ax.plot(kde_x, kde_y_scaled, color=color, linewidth=2.0, alpha=0.9)

    # Stat annotations

    median_val = np.median(arr)
    mean_val   = np.mean(arr)

    stddev_val = np.std(arr) 

    print(f" - For Field {field}, mean: {median_val:.4f}, std dev: {stddev_val:.4f}")

    if 0:
        ax.axvline(median_val, color="white", linewidth=1.4, linestyle="--",
                   label=f"median = {median_val:.2f}", alpha=0.9)
        ax.axvline(mean_val,   color="yellow", linewidth=1.2, linestyle=":",
                   label=f"mean = {mean_val:.2f}", alpha=0.8)
    
        
    ax.set_xlabel(format_axis_label(field), fontsize=11, color="#cccccc")
    ax.set_ylabel("Job count", fontsize=11, color="#cccccc")
    title = f"{field}"
    if jobname_pattern:
        title += f"  [jobs matching '{jobname_pattern}']"
    
    ax.set_title(title, fontsize=13, fontweight="bold", color="white", pad=10)

    ax.legend(fontsize=9, framealpha=0.3, labelcolor="white")

    # Stat box in corner
    stats_text = (
        f"n = {len(arr)}\n"
        f"min = {arr.min():.2f}\n"
        f"max = {arr.max():.2f}\n"
        f"σ = {arr.std():.2f}"
    )
    ax.text(
        0.98, 0.97, stats_text,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=8.5,
        color="#cccccc",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#1e1e2e", alpha=0.7, edgecolor="#444"),
    )


def apply_dark_style(fig: plt.Figure, axes) -> None:
    BG = "#000000" #"#12121e"
    PANEL = "#000000" #"#1c1c2e"
    fig.patch.set_facecolor(BG)
    for ax in (axes if hasattr(axes, "__iter__") else [axes]):
        ax.set_facecolor(PANEL)
        ax.tick_params(colors="#aaaaaa", labelsize=9)
        for spine in ax.spines.values():
            spine.set_edgecolor("#333355")
        ax.xaxis.label.set_color("#cccccc")
        ax.yaxis.label.set_color("#cccccc")


def make_figures(
    data: dict[str, list[float]],
    bins: int,
    jobname_pattern: str, # | None,
    outdir: str ,#| None,
) -> None:
    fields_with_data = [f for f, v in data.items() if v]
    fields_empty     = [f for f, v in data.items() if not v]

    if fields_empty:
        print(f"WARNING: No numeric data found for field(s): {', '.join(fields_empty)}")

    if not fields_with_data:
        sys.exit("ERROR: No plottable data found. Check your filters and field names.")

    n = len(fields_with_data)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(6.5 * ncols, 4.8 * nrows),
        squeeze=False,
    )

    for idx, field in enumerate(fields_with_data):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        color = PALETTE[idx % len(PALETTE)]
        plot_histogram(field, data[field], bins, ax, color, jobname_pattern)

    # Hide unused subplot slots
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    apply_dark_style(fig, axes.flat)

    fig.suptitle(
        "SLURM Job Statistics",
        fontsize=16, fontweight="bold", color="white", y=1.01,
    )
    fig.tight_layout()

    if outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        out_path = Path(outdir) / "slurm_hist.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        print(f"Saved → {out_path}")
    else:
        plt.show()


# ── CLI ────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Histogram wrapper for sacct SLURM job statistics.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # ── sacct filters ──
    grp = p.add_argument_group("sacct filter options")
    grp.add_argument("--user",       metavar="USER",
                     help="Filter by SLURM username (comma-separated list ok)")
    grp.add_argument("--starttime",  metavar="DATETIME",
                     help="Start of time window, e.g. 2025-01-01 or 2025-01-01T08:00")
    grp.add_argument("--endtime",    metavar="DATETIME",
                     help="End of time window")
    grp.add_argument("--account",    metavar="ACCOUNT",
                     help="Filter by SLURM account")
    grp.add_argument("--partition",  metavar="PARTITION",
                     help="Filter by partition name")
    grp.add_argument("--state",      metavar="STATE",
                     help="Job state(s), e.g. CD,F,TO  (default: CD = completed)")
    grp.add_argument("--jobs",       metavar="JOBID", nargs="+",
                     help="Explicit job IDs to query")
    grp.add_argument("--jobname",    metavar="PATTERN",
                     help="Glob wildcard filter on job name, e.g. 'sim_run_*'")

    # ── plot options ──
    grp2 = p.add_argument_group("plot options")
    grp2.add_argument("--fields",    metavar="FIELD", nargs="+", required=True,
                      help="sacct format fields to histogram, e.g. Elapsed MaxRSS CPUTime")
    grp2.add_argument("--bins",      type=int, default=20,
                      help="Number of histogram bins (default: 20)")
    grp2.add_argument("--outdir",    metavar="DIR",
                      help="Save figure to DIR/slurm_hist.png instead of displaying")

    # ── dry-run ──
    p.add_argument("--dry-run", action="store_true",
                   help="Print the sacct command that would be run, then exit")

    return p


def check_scipy() -> bool:
    try:
        import scipy.stats  # noqa: F401
        return True
    except ImportError:
        return False


def main() -> None:
    parser = build_parser()
    args   = parser.parse_args()

    if not check_scipy():
        print(
            "INFO: scipy not found — KDE overlay will be skipped. "
            "Install with: pip install scipy"
        )
        # Monkey-patch so the import inside plot_histogram fails gracefully
        import builtins, unittest.mock as mock  # noqa
        _real_import = builtins.__import__
        def _safe_import(name, *a, **kw):
            if name == "scipy.stats":
                raise ImportError
            return _real_import(name, *a, **kw)
        builtins.__import__ = _safe_import

    fields = args.fields
    cmd    = build_sacct_cmd(args, fields)

    if args.dry_run:
        print("Would run:", " ".join(cmd))
        sys.exit(0)

    print(f"Running: {' '.join(cmd)}")
    rows = run_sacct(cmd)

    if len(rows) < 1:
        print("sacct query returned no results.")
        return 

    print(f"sacct returned {len(rows)} row(s) — parsing…")

    data = collect_data(rows, fields, args.jobname)

    totals = {f: len(v) for f, v in data.items()}
    for f, n in totals.items():
        print(f"  {f}: {n} valid value(s)")

    make_figures(data, args.bins, args.jobname, args.outdir)


if __name__ == "__main__":
    main()
