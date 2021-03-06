#!/usr/bin/env python

## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

"""Combine directory information from multiple restarts to create unified
dumps.visit and lag_data.visit files. Call as

    ./make_visit_file.py

to combine files in the current working directory or as

    ./make_visit_file.py /path/to/other/directory

to combine files in some other directory.
"""
from __future__ import print_function

import os
from os.path import join
import sys

def find_files(base_file_name, directory):
    files = []
    if not directory:
        print("Directory not supplied. Using CWD")
        list_dir = os.listdir(os.getcwd())
    else:
        list_dir = os.listdir(directory)
    for f in list_dir:
        if base_file_name in f:
            files.append(f)
    new_files = []

    for num in range(16,25):
        temp_files = []
        for f in files:
            if len(f) == num:
                temp_files.append(f)
        temp_files.sort()
        if temp_files:
            new_files.extend(temp_files)

    if not len(new_files) == len(files):
        print("All files were not written.")
        print("num files found:  ", len(files))
        print("num files written:", len(new_files))
        print("Aborting")
        sys.exit(1)
    return new_files
    
num_args = len(sys.argv)
base_path = ""
if (num_args > 1):
    base_path = sys.argv[1]
else:
    base_path = os.getcwd()
print("Creating files in path:", base_path)

## Write Cartesian files
files = find_files("visit_dump.", base_path)
with open(base_path + '/dumps.visit', 'w') as f_vel:
    for f in files:
        f_vel.write(f + "/summary.samrai")
        f_vel.write("\n")

## Write Lagrangian files
files = find_files("lag_data.cycle_", base_path)
with open(base_path + '/lag_data.visit', 'w') as f_lag:
    for f in files:
        num_str = f[-6:]
        f_lag.write(f + "/lag_data.cycle_" + num_str + ".summary.silo")
        f_lag.write("\n")
