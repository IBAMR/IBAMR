if HAVE_FORTRAN
AM_CPPFLAGS = -I$(top_builddir) -I$(top_builddir)/include -I$(top_srcdir)/include
LDADD=libqdmod.a libqd_f_main.a $(top_builddir)/src/libqd.a $(FCLIBS)

if UPCASE_MODULE
DDEXT=DDEXT.$(module_ext)
QDEXT=QDEXT.$(module_ext)
QDMOD=QDMODULE.$(module_ext)
DDMOD=DDMODULE.$(module_ext)
else
DDEXT=ddext.$(module_ext)
QDEXT=qdext.$(module_ext)
QDMOD=qdmodule.$(module_ext)
DDMOD=ddmodule.$(module_ext)
endif

lib_LIBRARIES = libqdmod.a libqd_f_main.a
libqdmod_a_SOURCES = ddext.f ddmod.f qdext.f qdmod.f f_dd.cpp f_qd.cpp
libqd_f_main_a_SOURCES = main.cpp
ddmod.o: ddext.o
qdmod.o: ddmod.o qdext.o
$(QDMOD): qdmod.o
$(DDMOD): ddmod.o
$(DDEXT): ddext.o
$(QDEXT): qdext.o

pkglib_DATA = $(QDMOD) $(QDEXT) $(DDMOD) $(DDEXT)

DEMO=quaderq$(EXEEXT) quadgsq2d$(EXEEXT) quadgsq$(EXEEXT) \
     quadtsq$(EXEEXT) quadtsq2d$(EXEEXT)
TIMER = dd_timer$(EXEEXT) qd_timer$(EXEEXT)
EXTRA_PROGRAMS=$(DEMO) $(TIMER)

quaderq_SOURCES = tquaderq.f
nodist_quaderq_SOURCES=second.f
quaderq_LINK=$(CXXLINK)

quadgsq2d_SOURCES = tquadgsq2d.f
nodist_quadgsq2d_SOURCES=second.f
quadgsq2d_LINK=$(CXXLINK)

quadgsq_SOURCES = tquadgsq.f
nodist_quadgsq_SOURCES=second.f
quadgsq_LINK=$(CXXLINK)

quadtsq_SOURCES = tquadtsq.f
nodist_quadtsq_SOURCES=second.f
quadtsq_LINK=$(CXXLINK)

quadtsq2d_SOURCES = tquadtsq2d.f
nodist_quadtsq2d_SOURCES=second.f
quadtsq2d_LINK=$(CXXLINK)

dd_timer_SOURCES = dd_timer.f
nodist_dd_timer_SOURCES = second.f
dd_timer_LINK=$(CXXLINK)

qd_timer_SOURCES = qd_timer.f
nodist_qd_timer_SOURCES = second.f
qd_timer_LINK=$(CXXLINK)

demo: $(DEMO)

CLEANFILES = $(EXTRA_PROGRAMS) *.mod

endif

MAINTAINERCLEANFILES = dd_timer.f qd_timer.f
EXTRA_DIST = dd_timer.f qd_timer.f zz_timer.f Makefile.sample

dd_timer.f: zz_timer.f
	sed -e 's/@@/dd/' zz_timer.f >dd_timer.f

qd_timer.f: zz_timer.f
	sed -e 's/@@/qd/' zz_timer.f >qd_timer.f

if HAVE_FORTRAN
.PHONY: demo
endif

