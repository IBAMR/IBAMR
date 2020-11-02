#!/bin/bash
## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

# Replacement for cmake -E create_symlink that is much faster.
#
# $1 is the source test directory
# $2 is the build test directory

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# We have to avoid using loops or other such things to work with weird file
# names (e.g., paths that include spaces)
find "$INPUT_DIR" \( -name '*.input' -o -name '*.output' \) -exec ln -f -s {} "$OUTPUT_DIR" \;
