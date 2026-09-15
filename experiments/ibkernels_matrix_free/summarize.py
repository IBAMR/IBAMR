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
                 "fortran_over_cpp_median"))
for key, samples in sorted(groups.items()):
    median = statistics.median(samples)
    ratio = ""
    dimension, kernel, cells, markers, shuffled, operation = key
    if operation.startswith("cpp_loop"):
        peer = key[:-1] + (operation.replace("cpp_loop", "fortran_loop"),)
        if peer in groups:
            ratio = statistics.median(groups[peer]) / median
    elif operation.startswith("cpp_patch"):
        peer = key[:-1] + (operation.replace("cpp_patch", "LEInteractor"),)
        if peer in groups:
            ratio = statistics.median(groups[peer]) / median
    writer.writerow((*key, len(samples), median, min(samples), max(samples),
                     100 * statistics.stdev(samples) / statistics.mean(samples), ratio))
