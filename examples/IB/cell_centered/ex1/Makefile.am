## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C FeedbackForcer.C FeedbackForcer.h VelocityBcCoefs.C VelocityBcCoefs.h
EXTRA_DIST  = input2d.cylinder input3d.cylinder input3d.sphere
EXTRA_DIST += generate_cylinder2d.m generate_cylinder3d.m compute_St.m generate_sphere3d.C

EXTRA_DIST += cylinder2d_256.vertex cylinder2d_512.vertex cylinder2d_1024.vertex
EXTRA_DIST += cylinder2d_256.spring cylinder2d_512.spring cylinder2d_1024.spring
EXTRA_DIST += cylinder3d_128.vertex cylinder3d_256.vertex cylinder3d_512.vertex
EXTRA_DIST += sphere3d_64.vertex sphere3d_128.vertex sphere3d_256.vertex
EXTRA_DIST += sphere3d_64.spring sphere3d_128.spring sphere3d_256.spring

EXTRA_PROGRAMS =
if SAMRAI2D_ENABLED
EXTRA_PROGRAMS += main2d
endif
if SAMRAI3D_ENABLED
EXTRA_PROGRAMS += main3d
endif

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/lib/libIBAMR2d.a $(top_builddir)/lib/libIBAMR.a $(top_builddir)/ibtk/lib/libIBTK2d.a $(top_builddir)/ibtk/lib/libIBTK.a $(top_builddir)/ibtk/third_party/qd-2.3.7/fortran/libqdmod.a $(top_builddir)/ibtk/third_party/qd-2.3.7/src/libqd.a
main2d_SOURCES = $(SOURCES)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(top_builddir)/lib/libIBAMR3d.a $(top_builddir)/lib/libIBAMR.a $(top_builddir)/ibtk/lib/libIBTK3d.a $(top_builddir)/ibtk/lib/libIBTK.a $(top_builddir)/ibtk/third_party/qd-2.3.7/fortran/libqdmod.a $(top_builddir)/ibtk/third_party/qd-2.3.7/src/libqd.a
main3d_SOURCES = $(SOURCES)

examples: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input{2d.cylinder,3d.cylinder,3d.sphere} $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d_{256,512,1024}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/cylinder3d_{128,256,512}.vertex $(PWD) ; \
	  cp -f $(srcdir)/sphere3d_{64,128,256}.{vertex,spring} $(PWD) ; \
          cp -f $(srcdir)/compute_St.m $(PWD) ; \
	fi ;

clean-local:
	rm -f $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  rm -f $(builddir)/input{2d.cylinder,3d.cylinder,3d.sphere} ; \
	  rm -f $(builddir)/cylinder2d_{256,512,1024}.{vertex,spring} ; \
	  rm -f $(builddir)/cylinder3d_{128,256,512}.vertex ; \
	  rm -f $(builddir)/sphere3d_{64,128,256}.{vertex,spring} ; \
          rm -f $(builddir)/compute_St.m ; \
	fi ;
