## ---------------------------------------------------------------------
##
## Copyright (c) 2011 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

EXTRA_DIST = lagrangian_delta.f.m4 lagrangian_interaction2d.f.m4 lagrangian_interaction3d.f.m4
BUILT_SOURCES = lagrangian_delta.f lagrangian_interaction2d.f lagrangian_interaction3d.f
CLEANFILES = ${BUILT_SOURCES}

# we have to use some custom Make logic here since, long story short, autotools
# does not correctly generate Fortran file dependency information and, to get
# function inlining, we have to include lagrangian_delta.f in the other two
# files.
lagrangian_delta.f: lagrangian_delta.f.m4 Makefile
	$(M4) $(FM4FLAGS) $(AM_FM4FLAGS) -DTOP_SRCDIR=$(top_srcdir) -DCURRENT_SRCDIR=$(srcdir) \
    -DSAMRAI_FORTDIR=@SAMRAI_FORTDIR@ $< > $@

lagrangian_interaction2d.f: lagrangian_interaction2d.f.m4 lagrangian_delta.f Makefile
	$(M4) $(FM4FLAGS) $(AM_FM4FLAGS) -DTOP_SRCDIR=$(top_srcdir)  -DCURRENT_SRCDIR=$(srcdir) \
    -DSAMRAI_FORTDIR=@SAMRAI_FORTDIR@ $< > $@

lagrangian_interaction3d.f: lagrangian_interaction3d.f.m4 lagrangian_delta.f Makefile
	$(M4) $(FM4FLAGS) $(AM_FM4FLAGS) -DTOP_SRCDIR=$(top_srcdir)  -DCURRENT_SRCDIR=$(srcdir) \
    -DSAMRAI_FORTDIR=@SAMRAI_FORTDIR@ $< > $@
