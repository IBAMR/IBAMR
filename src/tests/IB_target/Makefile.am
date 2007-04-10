## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

## Dimension-dependent testers
SOURCES = main.C
EXTRA_DIST  = input2d.cylinder input3d.cylinder input3d.sphere
EXTRA_DIST += cylinder2d_256.vertex cylinder2d_512.vertex cylinder2d_1024.vertex
EXTRA_DIST += cylinder2d_256.spring cylinder2d_512.spring cylinder2d_1024.spring
EXTRA_DIST += cylinder3d_192.vertex
EXTRA_DIST += sphere3d_64.vertex sphere3d_128.vertex sphere3d_256.vertex
EXTRA_DIST += sphere3d_64.spring sphere3d_128.spring sphere3d_256.spring
EXTRA_DIST += generate_cylinder2d.m compute_St.m mesh_sphere.C

check_PROGRAMS = main2d main3d
check-local:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  cp -f $(srcdir)/input2d.cylinder input3d.cylinder input3d.sphere $(PWD) ; \
	  cp -f $(srcdir)/cylinder2d_{256,512,1024}.{vertex,spring} $(PWD) ; \
	  cp -f $(srcdir)/cylinder3d_192.vertex $(PWD) ; \
	  cp -f $(srcdir)/sphere3d_{64,128,256}.{vertex,spring} $(PWD) ; \
          cp -f $(srcdir)/compute_St.m $(PWD) ; \
	fi ;

main2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
main2d_LDADD = $(top_builddir)/src/libIBAMR2d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS2d.a $(top_builddir)/stools/src/libSTOOLS.a
main2d_SOURCES = $(SOURCES)

main3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
main3d_LDADD = $(top_builddir)/src/libIBAMR3d.a $(top_builddir)/src/libIBAMR.a $(top_builddir)/stools/src/libSTOOLS3d.a $(top_builddir)/stools/src/libSTOOLS.a
main3d_SOURCES = $(SOURCES)
