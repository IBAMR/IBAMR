#!/bin/bash

## ---------------------------------------------------------------------
##
## Copyright (c) 2018 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------


# Log file name
log_file=INS2d.log

# Number of lines to extract from the bottom of the log
num_lines=12

# Finest grid size
N=64

echo $PWD
echo "Extracting data..."

# Extra error data
tail -n $num_lines $log_file > temp_log
grep -oP ': \K.*' temp_log > temp_err
echo $N $(head -1 temp_err | tail -1) >> U_L1Error.txt
echo $N $(head -2 temp_err | tail -1) >> U_L2Error.txt
echo $N $(head -3 temp_err | tail -1) >> U_LInfError.txt
echo $N $(head -4 temp_err | tail -1) >> P_L1Error.txt
echo $N $(head -5 temp_err | tail -1) >> P_L2Error.txt
echo $N $(head -6 temp_err | tail -1) >> P_LInfError.txt

# Remove temporary files
rm temp_log
rm temp_err
