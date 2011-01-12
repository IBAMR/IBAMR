## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
doc:
	cd $(top_builddir)/ibtk/doc && $(MAKE) doc
	@DOXYGEN@

online_doc:
	cd $(top_builddir)/ibtk/doc && $(MAKE) online_doc
	@DOXYGEN@ online_Doxyfile
