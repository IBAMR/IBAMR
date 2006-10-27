## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = stools doc src
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = autogen.sh bin/includes doc/Doxyfile doc/libstdc++.tag include

## Automatically link header files in $(top_srcdir)/include and
## refresh emacs TAGS file.
##
## Note that these are for developer use only.
copy_headers:
	cd stools && $(MAKE) relink_headers
	cd $(top_srcdir) && bin/includes
relink_headers:
	cd stools && $(MAKE) relink_headers
	cd $(top_srcdir) && bin/includes --link
refresh_tags:
	cd stools && $(MAKE) refresh_tags
	cd $(top_srcdir) && find . -name "*.[ChIf]" -print | etags -
