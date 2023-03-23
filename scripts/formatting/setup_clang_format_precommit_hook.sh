#! /bin/bash
## ---------------------------------------------------------------------
##
## Copyright (c) 2016 - 2016 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

cp scripts/formatting/check_formatting.sh .git/hooks/pre-commit
chmod +x .git/hooks/pre-commit
