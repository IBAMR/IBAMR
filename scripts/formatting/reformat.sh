#! /bin/bash
## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2016 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

FILES="`find examples include ibtk/include ibtk/examples ibtk/src src tests -name *.cpp -o -name *.h`"
for file in $FILES; do
  echo $file
  clang-format -i -style=file $file
done
