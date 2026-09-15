#!/usr/bin/env python3
"""Summarize raw benchmark CSVs without discarding samples."""
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

groups = defaultdict(list)
for filename in sys.argv[1:]:
    with Path(filename).open() as stream:
        for row in csv.DictReader(stream):
            key = tuple(row[name] for name in
                        ("dimension", "kernel", "cells", "markers", "shuffled", "operation"))
            groups[key].append(float(row["seconds"]) / int(row["iterations"]))

writer = csv.writer(sys.stdout)
writer.writerow(("dimension", "kernel", "cells", "markers", "shuffled", "operation",
                 "samples", "median_seconds", "minimum_seconds", "maximum_seconds", "cv_percent",
                 "fortran_over_cpp_median", "expanded_over_cpp_median"))
for key, samples in sorted(groups.items()):
    median = statistics.median(samples)
    ratio = ""
    dimension, kernel, cells, markers, shuffled, operation = key
    expanded_ratio = ""
    if operation.startswith("cpp_"):
        scope, action = operation.split("_")[-2:]
        peer_operation = ("fortran_loop_" if scope == "loop" else "LEInteractor_") + action
        peer = key[:-1] + (peer_operation,)
        if peer in groups:
            ratio = statistics.median(groups[peer]) / median
        expanded = key[:-1] + ("cpp_expanded_" + scope + "_" + action,)
        if expanded in groups:
            expanded_ratio = statistics.median(groups[expanded]) / median
    writer.writerow((*key, len(samples), median, min(samples), max(samples),
                     100 * statistics.stdev(samples) / statistics.mean(samples), ratio, expanded_ratio))
