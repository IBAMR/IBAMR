#
# Written by Patrick Guio <patrick.guio@fys.uio.no>
#

arraydir = $(includedir)/blitz/array
generatedir = ../generate

genheaders = stencil-classes.cc 
#bops.cc uops.cc

array_HEADERS = asexpr.h asexpr.cc cartesian.h cgsolve.h complex.cc	\
convolve.cc convolve.h cycle.cc domain.h et.h expr.h expr.cc		\
fastiter.h funcs.h functorExpr.h geometry.h indirect.h interlace.cc	\
io.cc iter.h map.h methods.cc misc.cc multi.h newet-macros.h newet.h	\
ops.cc ops.h reduce.cc reduce.h resize.cc shape.h slice.h slicing.cc	\
stencil-et.h stencil-et-macros.h stencilops.h stencils.cc stencils.h	\
storage.h where.h zip.h $(genheaders)


clean-local:
	-rm -rf $(genheaders)
