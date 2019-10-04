#!/usr/bin/env python

## Filename: make_visit_file.py
##
## Copyright (c) 2002-2019, Boyce Griffith
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##    * Redistributions of source code must retain the above copyright notice,
##      this list of conditions and the following disclaimer.
##
##    * Redistributions in binary form must reproduce the above copyright
##      notice, this list of conditions and the following disclaimer in the
##      documentation and/or other materials provided with the distribution.
##
##    * Neither the name of The University of North Carolina nor the names of
##      its contributors may be used to endorse or promote products derived from
##      this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
## LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
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
