#!/bin/bash

# Log file name
log_file=INS2d.log

# Number of lines to extract from the bottom of the log
num_lines=18

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
echo $N $(head -7 temp_err | tail -1) >> RHO_L1Error.txt
echo $N $(head -8 temp_err | tail -1) >> RHO_L2Error.txt
echo $N $(head -9 temp_err | tail -1) >> RHO_LInfError.txt

# Remove temporary files
rm temp_log
rm temp_err
