#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

SUBDIRS = generate meta array

blitzdir = $(includedir)/blitz
generatedir = ./generate

genheaders = 
#matbops.h mathfunc.h matuops.h promote-old.h vecbops.cc vecuops.cc vecwhere.cc

blitz_HEADERS = array-impl.h array.cc array.h bench.cc bench.h	      \
benchext.cc benchext.h blitz.h bounds.h bzconfig.h bzdebug.h	      \
compiler.h constpointerstack.h etbase.h et-forward.h funcs.h	      \
globeval.cc indexexpr.h indexmap-forward.h levicivita.h limits-hack.h \
listinit.h memblock.cc memblock.h minmax.h numinquire.h numtrait.h    \
ops.h prettyprint.h promote.h range.h range.cc ranks.h reduce.h	      \
shapecheck.h simdtypes.h tau.h timer.h tinymat2.h tinymat2.cc	      \
tinymat2io.cc tinyvec2.cc tinyvec2.h tinyvec2io.cc tm2fastiter.h      \
tmevaluate.h tv2fastiter.h tvevaluate.h traversal.cc traversal.h      \
tuning.h tvcross.h tvecglobs.h update.h vector2.h wrap-climits.h      \
$(genheaders)

EXTRA_HEADERS = apple/bzconfig.h intel/bzconfig.h ibm/bzconfig.h \
compaq/bzconfig.h hp/bzconfig.h sgi/bzconfig.h gnu/bzconfig.h \
pgi/bzconfig.h pathscale/bzconfig.h kai/bzconfig.h fujitsu/bzconfig.h

nobase_blitz_HEADERS = $(COMPILER_SPECIFIC_HEADER)

DISTCLEANFILES = apple/bzconfig.h intel/bzconfig.h ibm/bzconfig.h \
compaq/bzconfig.h hp/bzconfig.h sgi/bzconfig.h gnu/bzconfig.h \
pgi/bzconfig.h pathscale/bzconfig.h kai/bzconfig.h fujitsu/bzconfig.h

clean-local:
	-rm -rf config.h $(genheaders)
