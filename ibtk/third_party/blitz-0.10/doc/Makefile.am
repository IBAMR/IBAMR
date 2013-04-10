#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

SUBDIRS = examples stencils doxygen

# docdir directory to install the doc (pdf ps html) 
docdir = $(datadir)/doc/$(PACKAGE)-$(VERSION)

DVIPS=export TEXPICTS=$(srcdir); dvips
AM_MAKEINFOFLAGS = --no-split -I$(srcdir)
AM_MAKEINFOHTMLFLAGS = -I$(srcdir)


STENCILS = stencils/backward11.texi stencils/backward12.texi \
stencils/backward21.texi stencils/backward22.texi stencils/backward31.texi \
stencils/backward32.texi stencils/backward41.texi stencils/backward42.texi \
stencils/central12.texi stencils/central14.texi stencils/central22.texi \
stencils/central24.texi stencils/central32.texi stencils/central34.texi \
stencils/central42.texi stencils/central44.texi stencils/forward11.texi \
stencils/forward12.texi stencils/forward21.texi stencils/forward22.texi \
stencils/forward31.texi stencils/forward32.texi stencils/forward41.texi \
stencils/forward42.texi stencils/Laplacian2D4.texi stencils/Laplacian2D.texi

EXAMPLES = examples/cast.texi examples/debug.texi examples/dump.texi \
examples/fixed-point.texi examples/fixed.texi examples/io.texi \
examples/outer.texi examples/output.texi examples/range.texi \
examples/simple.texi examples/slicing.texi examples/storage.texi \
examples/strideslice.texi examples/xor.texi

OUTPUT = examples/cast.out examples/debug.out examples/dump.out \
examples/fixed.out examples/io.out examples/outer.out examples/output.out \
examples/range.out examples/simple.out examples/slicing.out \
examples/storage.out examples/strideslice.out examples/xor.out

# if program makeinfo available 
# setup _TEXINFOS target
if MAKEINFO_DOC

INFO_SRC = blitz.texi
BLITZ_SRC = about.texi arrays-ctors.texi arrays-debug.texi \
arrays-expr.texi arrays-globals.texi arrays-indirect.texi arrays-intro.texi \
arrays-io.texi arrays-members.texi arrays-multi.texi arrays-slicing.texi \
arrays-stencils.texi arrays-storage.texi arrays-types.texi \
arrays-usertype.texi compiling.texi constants.texi download.texi \
faq.texi help.texi install.texi legal.texi numinquire.texi parallel.texi \
platforms.texi random.texi tau.texi tinymatrix.texi tinyvector.texi \
tuning.texi copyright.texi $(STENCILS) $(EXAMPLES) $(OUTPUT)

else

INFO_SRC = 
BLITZ_SRC = 

endif

info_TEXINFOS = $(INFO_SRC) 
blitz_TEXINFOS = $(BLITZ_SRC)

EXTRA_DIST = blitz.info blitz.dvi blitz.ps blitz.pdf blitz.gif blitztiny.jpg \
	indirect.eps indirect.fig indirect.pdf indirect.gif indirect.txt \
	sinsoid.eps sinsoid.pdf sinsoid.gif sinsoid.txt \
	slice.eps slice.fig slice.pdf slice.gif slice.txt \
	strideslice.eps strideslice.fig strideslice.pdf strideslice.gif strideslice.txt \
	tensor1.eps tensor1.fig tensor1.pdf tensor1.gif tensor1.txt 



# if program texi2html available 
# use it to generate html documentation
if TEXI2HTML_DOC

TEXI2HTML_FLAGS = --subdir=blitz.html --split=chapter

blitz.html: $(info_TEXINFOS) $(blitz_TEXINFOS)
	$(TEXI2HTML) $(TEXI2HTML_FLAGS) $(srcdir)/blitz.texi

# Install also the images into the directory html 
install-html-local:
	test -z "$(htmldir)" || $(MKDIR_P) "$(DESTDIR)$(htmldir)"
	$(INSTALL_DATA) *.gif *.jpg $(DESTDIR)$(htmldir)

.PHONY: html

else

# Install also the images into the directory containing the html files
install-html-local:
	$(INSTALL_DATA) *.gif *.jpg blitz.html

endif

# extra dependencies 
stencils/%.texi:
	(cd stencils; $(MAKE) $(AM_MAKEFLAGS) stencils)

examples/%.texi:
	(cd examples; $(MAKE) $(AM_MAKEFLAGS) -i texi)

examples/%.out:
	(cd examples; $(MAKE) $(AM_MAKEFLAGS) -i out)


# Files from texi2dvi that should be removed, but which Automake does
# not know.
CLEANFILES = blitz.tmp

clean-local:
	-rm -rf blitz.html 

