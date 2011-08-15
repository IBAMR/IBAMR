## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

lib_LIBRARIES = libIBAMR.a
if SAMRAI2D_ENABLED
lib_LIBRARIES += libIBAMR2d.a
endif
if SAMRAI3D_ENABLED
lib_LIBRARIES += libIBAMR3d.a
endif

all-local: $(lib_LIBRARIES)
	for file in "lib*.a" ; do \
	  ln -fs $(PWD)/$$file $(PWD)/../src; \
	done

clean-local:
	rm -f $(PWD)/../src/libIBAMR.a
	rm -f $(PWD)/../src/libIBAMR2d.a
	rm -f $(PWD)/../src/libIBAMR3d.a

## Dimension-independent library
libIBAMR_a_SOURCES = \
$(top_builddir)/src/advect/fortran/godunov_helpers.f \
$(top_builddir)/src/fortran/sqrtm.f

pkg_include_HEADERS = \
$(top_srcdir)/src/utilities/ibamr_enums.h \
$(top_srcdir)/src/utilities/ibamr_utilities.h \
$(top_srcdir)/src/utilities/namespaces.h

## Dimension-dependent libraries

###
###$(top_srcdir)/src/IB/IBHierarchyIntegrator.C \
###

DIM_INDEPENDENT_SOURCES = \
$(top_srcdir)/src/IB/IBAnchorPointSpec.C \
$(top_srcdir)/src/IB/IBAnchorPointSpecFactory.C \
$(top_srcdir)/src/IB/IBBeamForceSpec.C \
$(top_srcdir)/src/IB/IBBeamForceSpecFactory.C \
$(top_srcdir)/src/IB/IBEulerianForceFunction.C \
$(top_srcdir)/src/IB/IBEulerianSourceFunction.C \
$(top_srcdir)/src/IB/IBHDF5Initializer.C \
$(top_srcdir)/src/IB/IBInstrumentPanel.C \
$(top_srcdir)/src/IB/IBInstrumentationSpec.C \
$(top_srcdir)/src/IB/IBInstrumentationSpecFactory.C \
$(top_srcdir)/src/IB/IBKirchhoffRodForceGen.C \
$(top_srcdir)/src/IB/IBLagrangianForceStrategy.C \
$(top_srcdir)/src/IB/IBLagrangianForceStrategySet.C \
$(top_srcdir)/src/IB/IBLagrangianSourceStrategy.C \
$(top_srcdir)/src/IB/IBPostProcessStrategy.C \
$(top_srcdir)/src/IB/IBRodForceSpec.C \
$(top_srcdir)/src/IB/IBRodForceSpecFactory.C \
$(top_srcdir)/src/IB/IBSourceSpec.C \
$(top_srcdir)/src/IB/IBSourceSpecFactory.C \
$(top_srcdir)/src/IB/IBSpringForceSpec.C \
$(top_srcdir)/src/IB/IBSpringForceSpecFactory.C \
$(top_srcdir)/src/IB/IBStandardForceGen.C \
$(top_srcdir)/src/IB/IBStandardSourceGen.C \
$(top_srcdir)/src/IB/IBStandardInitializer.C \
$(top_srcdir)/src/IB/IBTargetPointForceSpec.C \
$(top_srcdir)/src/IB/IBTargetPointForceSpecFactory.C \
$(top_srcdir)/src/adv_diff/AdvDiffHierarchyIntegrator.C \
$(top_srcdir)/src/adv_diff/AdvDiffHypPatchOps.C \
$(top_srcdir)/src/advect/AdvectHypPatchOps.C \
$(top_srcdir)/src/advect/GodunovAdvector.C \
$(top_srcdir)/src/base_hier_integrator/HierarchyIntegrator.C \
$(top_srcdir)/src/navier_stokes/INSHierarchyIntegrator.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredHierarchyIntegrator.C \
$(top_srcdir)/src/navier_stokes/INSProblemCoefs.C

###if LIBMESH_ENABLED
###DIM_INDEPENDENT_SOURCES += \
###$(top_srcdir)/src/IB/IBFEHierarchyIntegrator.C
###endif

###
###$(top_srcdir)/src/IB/IBHierarchyIntegrator.h \
###

pkg_include_HEADERS += \
$(top_srcdir)/src/IB/IBAnchorPointSpec.I \
$(top_srcdir)/src/IB/IBAnchorPointSpec.h \
$(top_srcdir)/src/IB/IBAnchorPointSpecFactory.h \
$(top_srcdir)/src/IB/IBBeamForceSpec.I \
$(top_srcdir)/src/IB/IBBeamForceSpec.h \
$(top_srcdir)/src/IB/IBBeamForceSpecFactory.h \
$(top_srcdir)/src/IB/IBEulerianForceFunction.h \
$(top_srcdir)/src/IB/IBEulerianSourceFunction.h \
$(top_srcdir)/src/IB/IBHDF5Initializer.h \
$(top_srcdir)/src/IB/IBInstrumentPanel.h \
$(top_srcdir)/src/IB/IBInstrumentationSpec.I \
$(top_srcdir)/src/IB/IBInstrumentationSpec.h \
$(top_srcdir)/src/IB/IBInstrumentationSpecFactory.h \
$(top_srcdir)/src/IB/IBKirchhoffRodForceGen.h \
$(top_srcdir)/src/IB/IBLagrangianForceStrategy.h \
$(top_srcdir)/src/IB/IBLagrangianForceStrategySet.h \
$(top_srcdir)/src/IB/IBLagrangianSourceStrategy.h \
$(top_srcdir)/src/IB/IBPostProcessStrategy.h \
$(top_srcdir)/src/IB/IBRodForceSpec.I \
$(top_srcdir)/src/IB/IBRodForceSpec.h \
$(top_srcdir)/src/IB/IBRodForceSpecFactory.h \
$(top_srcdir)/src/IB/IBSourceSpec.I \
$(top_srcdir)/src/IB/IBSourceSpec.h \
$(top_srcdir)/src/IB/IBSourceSpecFactory.h \
$(top_srcdir)/src/IB/IBSpringForceFunctions.h \
$(top_srcdir)/src/IB/IBSpringForceSpec.I \
$(top_srcdir)/src/IB/IBSpringForceSpec.h \
$(top_srcdir)/src/IB/IBSpringForceSpecFactory.h \
$(top_srcdir)/src/IB/IBStandardForceGen.h \
$(top_srcdir)/src/IB/IBStandardSourceGen.h \
$(top_srcdir)/src/IB/IBStandardInitializer.h \
$(top_srcdir)/src/IB/IBTargetPointForceSpec.I \
$(top_srcdir)/src/IB/IBTargetPointForceSpec.h \
$(top_srcdir)/src/IB/IBTargetPointForceSpecFactory.h \
$(top_srcdir)/src/adv_diff/AdvDiffHierarchyIntegrator.h \
$(top_srcdir)/src/adv_diff/AdvDiffHypPatchOps.h \
$(top_srcdir)/src/advect/AdvectHypPatchOps.h \
$(top_srcdir)/src/advect/GodunovAdvector.h \
$(top_srcdir)/src/base_hier_integrator/HierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSHierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredHierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSProblemCoefs.h

###if LIBMESH_ENABLED
###pkg_include_HEADERS += \
###$(top_srcdir)/src/IB/IBFEHierarchyIntegrator.h
###endif

libIBAMR2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR2d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/adv_diff/fortran/adv_diff_consdiff2d.f \
$(top_builddir)/src/advect/fortran/advect_detect2d.f \
$(top_builddir)/src/advect/fortran/advect_diff2d.f \
$(top_builddir)/src/advect/fortran/advect_stable2d.f \
$(top_builddir)/src/advect/fortran/godunov_predictors2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_bdryop2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_divsource2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stabledt2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_derivatives2d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_helpers2d.f

libIBAMR3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBAMR3d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/adv_diff/fortran/adv_diff_consdiff3d.f \
$(top_builddir)/src/advect/fortran/advect_detect3d.f \
$(top_builddir)/src/advect/fortran/advect_diff3d.f \
$(top_builddir)/src/advect/fortran/advect_stable3d.f \
$(top_builddir)/src/advect/fortran/godunov_predictors3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_bdryop3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_divsource3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stabledt3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_derivatives3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_helpers3d.f
