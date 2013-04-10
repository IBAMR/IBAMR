#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

EXTRA_DIST = fixed-class.cpp fixed-point.h makefile.example io.data

AM_CPPFLAGS = -I$(srcdir) -I$(top_srcdir) -I$(top_builddir)
AM_CXXFLAGS = @CXX_DEBUG_FLAGS@
LDADD = -L$(top_builddir)/lib -lblitz

EXAMPLES = simple storage outer strideslice slicing debug output io dump xor \
cast range fixed

EXTRA_PROGRAMS =  $(EXAMPLES)

simple_SOURCES = simple.cpp
storage_SOURCES = storage.cpp
outer_SOURCES = outer.cpp
strideslice_SOURCES = strideslice.cpp
slicing_SOURCES = slicing.cpp
debug_SOURCES = debug.cpp
output_SOURCES = output.cpp
io_SOURCES = io.cpp
dump_SOURCES = dump.cpp
xor_SOURCES = xor.cpp
cast_SOURCES = cast.cpp
range_SOURCES = range.cpp
fixed_SOURCES = fixed.cpp

%.texi : 
	sed -e 's/{/@{/' -e 's/}/@}/' $< | fold -s -w 80 > $@

%.out : 
	$(MAKE) $(AM_MAKEFLAGS) $(patsubst %.out,%$(EXEEXT),$@)
	./$(patsubst %.out,%$(EXEEXT),$@) 2>&1 | fold -s -w 80 > $@

simple.texi: simple.cpp
storage.texi: storage.cpp
outer.texi: outer.cpp
strideslice.texi: strideslice.cpp
slicing.texi: slicing.cpp
debug.texi: debug.cpp
output.texi: output.cpp
io.texi: io.cpp
dump.texi: dump.cpp
xor.texi: xor.cpp
cast.texi: cast.cpp
range.texi: range.cpp
fixed.texi: fixed.cpp
fixed-point.texi: fixed-point.h

simple.out: simple.cpp
storage.out: storage.cpp
outer.out: outer.cpp
strideslice.out: strideslice.cpp
slicing.out: slicing.cpp
debug.out: debug.cpp
output.out: output.cpp
io.out: io.cpp
dump.out: dump.cpp
xor.out: xor.cpp
cast.out: cast.cpp
range.out: range.cpp
fixed.out: fixed.cpp

texi: simple.texi storage.texi outer.texi strideslice.texi slicing.texi \
debug.texi output.texi io.texi dump.texi xor.texi cast.texi range.texi \
fixed.texi fixed-point.texi

out: simple.out storage.out outer.out strideslice.out slicing.out \
debug.out output.out io.out dump.out xor.out cast.out range.out fixed.out


clean-local:
	-rm -rf *.ii *.ti cxx_repository Template.dir ii_files ti_files $(EXTRA_PROGRAMS) *.texi *.out core.[0-9]*


