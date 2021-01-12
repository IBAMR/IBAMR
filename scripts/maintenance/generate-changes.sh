#!/bin/bash
## ---------------------------------------------------------------------
##
## Copyright (c) 2019 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

if test ! -d ibtk -o ! -d include -o ! -d examples ; then
  echo "*** This script must be run from the top-level directory of IBAMR."
  exit 1
fi

cd doc/news/changes/

echo "<h2>Incompatibilities</h2>"
echo "<ol>"
for FILE in $(ls incompatibilities/2*)
do
    echo "<li>"
    cat $FILE
    echo "</li>"
done
echo "</ol>"

echo ""
echo "<h2>Major Changes</h2>"
echo "<ol>"
for FILE in $(ls major/2*)
do
    echo "<li>"
    cat $FILE
    echo "</li>"
done
echo "</ol>"

echo ""
echo "<h2>Minor Changes</h2>"
echo "<ol>"
for FILE in $(ls minor/2*)
do
    echo "<li>"
    cat $FILE
    echo "</li>"
done
echo "</ol>"
