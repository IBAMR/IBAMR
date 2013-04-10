#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

noinst_HEADERS = arroperands.h arroptuple.h bzfstream.h operands.h operands2.h optuple.h optuple2.h

EXTRA_DIST = genstencils.py

EXTRA_PROGRAMS = 

$(top_srcdir)/blitz/array/stencil-classes.cc: genstencils.py
	python $< $@

generate-headers: $(top_srcdir)/blitz/array/stencil-classes.cc

all-local: $(EXTRA_PROGRAMS)
	$(MAKE) $(AM_MAKEFLAGS) generate-headers

clean-local:
	-rm -rf cxx_repository ti_files ii_files $(EXTRA_PROGRAMS)


%.h: 
	./$(patsubst %.cpp,%$(EXEEXT),$<)

%.cc: 
	./$(patsubst %.cpp,%$(EXEEXT),$<)


