#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

EXTRA_DIST = dump-stencil.cpp

AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_builddir)
LDADD = -L$(top_builddir)/lib -lblitz
COMPILE = $(CXX) $(DEFS) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)

LTCXXCOMPILE = $(LIBTOOL) --mode=link $(CXX) $(DEFS) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS)

%2D.texi :
	$(LTCXXCOMPILE) -Dstenciltoapply=$(patsubst %.texi,%_stencilop,$@) \
		-Dstencilargs=A -o dump-stencil$(EXEEXT) $< $(LDADD)
	./dump-stencil$(EXEEXT) > $@  2>&1
	-rm -f dump-stencil$(EXEEXT)

%2D4.texi :
	$(LTCXXCOMPILE) -Dstenciltoapply=$(patsubst %.texi,%_stencilop,$@) \
		-Dstencilargs=A -o dump-stencil$(EXEEXT) $< $(LDADD)
	./dump-stencil$(EXEEXT) > $@  2>&1
	-rm -f dump-stencil$(EXEEXT)

%.texi :
	$(LTCXXCOMPILE) -Dstenciltoapply=$(patsubst %.texi,%_stencilop,$@) \
		-Dstencilargs=A,1 -o dump-stencil$(EXEEXT) $< $(LDADD)
	./dump-stencil$(EXEEXT) > $@  2>&1
	-rm -f dump-stencil$(EXEEXT)

central12.texi: dump-stencil.cpp
central22.texi: dump-stencil.cpp
central32.texi: dump-stencil.cpp
central42.texi: dump-stencil.cpp
central14.texi: dump-stencil.cpp
central24.texi: dump-stencil.cpp
central34.texi: dump-stencil.cpp
central44.texi: dump-stencil.cpp

forward11.texi: dump-stencil.cpp
forward21.texi: dump-stencil.cpp
forward31.texi: dump-stencil.cpp
forward41.texi: dump-stencil.cpp
forward12.texi: dump-stencil.cpp
forward22.texi: dump-stencil.cpp
forward32.texi: dump-stencil.cpp
forward42.texi: dump-stencil.cpp

backward11.texi: dump-stencil.cpp
backward21.texi: dump-stencil.cpp
backward31.texi: dump-stencil.cpp
backward41.texi: dump-stencil.cpp
backward12.texi: dump-stencil.cpp
backward22.texi: dump-stencil.cpp
backward32.texi: dump-stencil.cpp
backward42.texi: dump-stencil.cpp

Laplacian2D.texi: dump-stencil.cpp
Laplacian2D4.texi: dump-stencil.cpp

stencils: central12.texi central22.texi central32.texi central42.texi \
central14.texi central24.texi central34.texi central44.texi forward11.texi \
forward21.texi forward31.texi forward41.texi forward12.texi forward22.texi \
forward32.texi forward42.texi backward11.texi backward21.texi \
backward31.texi backward41.texi backward12.texi backward22.texi \
backward32.texi backward42.texi Laplacian2D.texi Laplacian2D4.texi


clean-local:
	-rm -rf *.ii *.ti cxx_repository Template.dir ii_files ti_files dump-stencil *.texi

