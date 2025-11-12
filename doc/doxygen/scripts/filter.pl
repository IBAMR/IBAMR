#!/usr/bin/perl

## ---------------------------------------------------------------------
##
## Copyright (c) 2025 - 2025 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

# This script is a subset of the filters used by the equivalent script in
# deal.II.

while (<>)
{
    # Convert SAMRAI-style \f$ to @f$. Explanation:
    #
    # (?<!)       : negative lookbehind assertion
    # (?<!\\)     : ensure that the previous character is not a '\'
    # (?!)        : negative lookahead assertion
    # (?!\$)      : ensure that the next character is not a '$'
    s#(?<!\\)\\f\$(?!\$)#\@f\$#g;

    # Fix links. This assumes that we don't have multiline markdown links to
    # modules
    s#.Style Guide..style-guide.md.#\@ref styleguide #g;
    s#.Project Architecture..architecture.md.#\@ref architecture #g;

    print;
}
