#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

ACLOCAL_AMFLAGS  = -I m4

if FORTRAN_ENABLED
FORTRAN_DIRS = benchmarks
endif

SUBDIRS = bin blitz random lib testsuite \
	examples $(FORTRAN_DIRS) doc

EXTRA_DIST = COPYRIGHT LEGAL LICENSE README.binutils README-VS.NET.txt \
	Blitz-VS.NET.zip Blitz-VS2005.NET.zip Blitz-VS2010.zip \
	ChangeLog.0 ChangeLog.1 ChangeLog.2

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = blitz.pc

DISTCLEANFILES = blitz.pc blitz.pc.in blitz-uninstalled.pc

lib: blitz-library

blitz-library: 
	@(cd blitz && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	@(cd lib && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

check-testsuite check-examples check-benchmarks: lib
	@(cd $(subst check-,,$@) && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;

clean-local:
	-rm -rf *.ii *.ti cxx_repository Template.dir ii_files ti_files

