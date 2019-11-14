include $(top_srcdir)/config/Make-rules

EXTRA_PROGRAMS = cf_relaxation_op_01_2d cf_relaxation_op_01_3d cf_forcing_op_01_2d cf_forcing_op_01_3d cf_four_roll_mill

cf_relaxation_op_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
cf_relaxation_op_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
cf_relaxation_op_01_2d_SOURCES = cf_relaxation_op_01.cpp

cf_relaxation_op_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
cf_relaxation_op_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
cf_relaxation_op_01_3d_SOURCES = cf_relaxation_op_01.cpp

cf_forcing_op_01_2d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
cf_forcing_op_01_2d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
cf_forcing_op_01_2d_SOURCES = cf_forcing_op_01.cpp

cf_forcing_op_01_3d_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
cf_forcing_op_01_3d_LDADD = $(IBAMR_LDFLAGS) $(IBAMR3d_LIBS) $(IBAMR_LIBS)
cf_forcing_op_01_3d_SOURCES = cf_forcing_op_01.cpp

cf_four_roll_mill_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
cf_four_roll_mill_LDADD = $(IBAMR_LDFLAGS) $(IBAMR2d_LIBS) $(IBAMR_LIBS)
cf_four_roll_mill_SOURCES = cf_four_roll_mill.cpp

tests: $(EXTRA_PROGRAMS)
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
	  ln -f -s $(srcdir)/*input $(PWD) ; \
	  ln -f -s $(srcdir)/*output $(PWD) ; \
	fi ;
.PHONY: tests
