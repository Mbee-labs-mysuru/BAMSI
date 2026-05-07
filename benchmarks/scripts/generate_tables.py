#!/usr/bin/env python3
"""BAMSI Benchmark Table Generator — Execution Plan §5.3.4

Reads raw CSV results from benchmarks/results/v1.0.0/ and produces
publication-ready markdown tables for Paper 1 inclusion.

Usage: python3 benchmarks/scripts/generate_tables.py [--dir results/v1.0.0/]
"""

import csv
import glob
import os
import sys
from collections import defaultdict

def load_csvs(results_dir):
    """Load all CSV benchmark results."""
    rows = []
    for path in sorted(glob.glob(os.path.join(results_dir, "*.csv"))):
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows.append(row)
    return rows

def aggregate(rows, group_keys, value_key):
    """Aggregate values by group keys, computing mean and stddev."""
    groups = defaultdict(list)
    for row in rows:
        key = tuple(row[k] for k in group_keys)
        try:
            groups[key].append(float(row[value_key]))
        except (ValueError, KeyError):
            pass

    result = {}
    for key, values in groups.items():
        n = len(values)
        mean = sum(values) / n if n else 0
        if n > 1:
            variance = sum((v - mean) ** 2 for v in values) / (n - 1)
            stddev = variance ** 0.5
        else:
            stddev = 0
        result[key] = {"mean": mean, "stddev": stddev, "n": n}
    return result

def generate_compression_table(rows):
    """Generate compression ratio table."""
    build_rows = [r for r in rows if r["operation"] == "build"]
    agg = aggregate(build_rows, ["dataset", "tool"], "wall_time_s")
    size_agg = aggregate(build_rows, ["dataset", "tool"], "file_size_bytes")

    if not agg:
        return "No compression data available.\n"

    lines = ["## Compression Results\n"]
    lines.append("| Dataset | Tool | Build Time (s) | Compressed Size (bytes) |")
    lines.append("|---------|------|---------------|------------------------|")

    for (dataset, tool), stats in sorted(agg.items()):
        size_stats = size_agg.get((dataset, tool), {"mean": 0, "stddev": 0})
        lines.append(f"| {dataset} | {tool} | "
                     f"{stats['mean']:.3f} ± {stats['stddev']:.3f} | "
                     f"{size_stats['mean']:.0f} |")

    return "\n".join(lines) + "\n"

def generate_query_table(rows):
    """Generate query latency table."""
    query_rows = [r for r in rows if r["operation"].startswith("count_")]
    agg = aggregate(query_rows, ["dataset", "tool", "operation"], "wall_time_s")

    if not agg:
        return "No query data available.\n"

    lines = ["\n## Query Latency Results\n"]
    lines.append("| Dataset | Tool | Pattern | Latency (s) |")
    lines.append("|---------|------|---------|-------------|")

    for (dataset, tool, op), stats in sorted(agg.items()):
        pattern = op.replace("count_", "")
        lines.append(f"| {dataset} | {tool} | {pattern} | "
                     f"{stats['mean']:.6f} ± {stats['stddev']:.6f} |")

    return "\n".join(lines) + "\n"

def generate_verify_reconstruct_table(rows):
    """Generate verify/reconstruct table."""
    vr_rows = [r for r in rows if r["operation"] in ("verify", "reconstruct")]
    agg = aggregate(vr_rows, ["dataset", "tool", "operation"], "wall_time_s")

    if not agg:
        return "No verify/reconstruct data available.\n"

    lines = ["\n## Verify & Reconstruct Results\n"]
    lines.append("| Dataset | Tool | Operation | Time (s) |")
    lines.append("|---------|------|-----------|----------|")

    for (dataset, tool, op), stats in sorted(agg.items()):
        lines.append(f"| {dataset} | {tool} | {op} | "
                     f"{stats['mean']:.3f} ± {stats['stddev']:.3f} |")

    return "\n".join(lines) + "\n"

def main():
    results_dir = "benchmarks/results/v1.0.0"
    if len(sys.argv) > 2 and sys.argv[1] == "--dir":
        results_dir = sys.argv[2]

    rows = load_csvs(results_dir)
    if not rows:
        print(f"No CSV files found in {results_dir}")
        sys.exit(1)

    output = "# BAMSI Benchmark Results — v1.0.0\n\n"
    output += f"Generated from {len(rows)} data points.\n\n"
    output += generate_compression_table(rows)
    output += generate_query_table(rows)
    output += generate_verify_reconstruct_table(rows)

    output_file = os.path.join(results_dir, "benchmark_tables.md")
    with open(output_file, "w") as f:
        f.write(output)

    print(output)
    print(f"\nTables written to: {output_file}")

if __name__ == "__main__":
    main()
