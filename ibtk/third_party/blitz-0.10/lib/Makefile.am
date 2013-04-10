#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

EXTRA_DIST = readme.txt

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir) $(BOOST_CPPFLAGS)

lib_LIBRARIES = libblitz.a
libblitz_a_SOURCES = $(top_srcdir)/src/globals.cpp
###libblitz_a_LIBADD = $(LIBS) $(CXX_LIBS)


clean-local:
	-rm -rf *.ii *.ti cxx_repository Template.dir ii_files ti_files
