## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = src
ACLOCAL_AMFLAGS = -I m4

## Automatically link header files in $(top_srcdir)/include/ibamr
BUILT_SOURCES = relink_headers
relink_headers:
	cd $(top_srcdir) && bin/includes --link
EXTRA_DIST = autogen.sh bin/includes doc/Doxyfile
dist-hook:
	mkdir $(distdir)/include
	mkdir $(distdir)/include/@PACKAGE@
