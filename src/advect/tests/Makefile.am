## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

EXTRA_DIST = input1d input2d input3d

## Dimension-dependent testers
SOURCES = main.C QInit.C QInit.h USet.C USet.h

check_PROGRAMS = main1d main2d main3d

main1d_SOURCES = $(SOURCES)
main1d_LDADD = ../../libIBAMR1d.la ../../libIBAMR.la
main1d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=1

main2d_SOURCES = $(SOURCES)
main2d_LDADD = ../../libIBAMR2d.la ../../libIBAMR.la
main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2

main3d_SOURCES = $(SOURCES)
main3d_LDADD = ../../libIBAMR3d.la ../../libIBAMR.la
main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
