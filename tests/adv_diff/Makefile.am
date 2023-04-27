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

include $(top_srcdir)/config/Make-rules


EXTRA_PROGRAMS = adv_diff_01_3d adv_diff_02_2d adv_diff_02_3d adv_diff_03_2d adv_diff_convec_opers_2d adv_diff_convec_opers_3d adv_diff_regridding_2d bp_adv_diff_01_2d bp_adv_diff_02_2d

adv_diff_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
adv_diff_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
adv_diff_01_3d_SOURCES = adv_diff_01.cpp

adv_diff_02_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
adv_diff_02_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
adv_diff_02_2d_SOURCES = adv_diff_02.cpp

adv_diff_02_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
adv_diff_02_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
adv_diff_02_3d_SOURCES = adv_diff_02.cpp

adv_diff_03_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
adv_diff_03_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
adv_diff_03_2d_SOURCES = adv_diff_03.cpp

adv_diff_convec_opers_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
adv_diff_convec_opers_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
adv_diff_convec_opers_2d_SOURCES = adv_diff_convec_opers.cpp

adv_diff_convec_opers_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
adv_diff_convec_opers_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
adv_diff_convec_opers_3d_SOURCES = adv_diff_convec_opers.cpp

adv_diff_regridding_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
adv_diff_regridding_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
adv_diff_regridding_2d_SOURCES = adv_diff_regridding.cpp

bp_adv_diff_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
bp_adv_diff_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
bp_adv_diff_01_2d_SOURCES = bp_adv_diff_01.cpp

bp_adv_diff_02_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
bp_adv_diff_02_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
bp_adv_diff_02_2d_SOURCES = bp_adv_diff_02.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
