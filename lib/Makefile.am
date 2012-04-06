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
$(top_builddir)/src/fortran/sqrtm.f \
$(top_srcdir)/src/navier_stokes/INSProblemCoefs.C

if LIBMESH_ENABLED
libIBAMR_a_SOURCES += \
$(top_srcdir)/src/solid_mechanics/ExplicitFEMechanicsSolver.C
endif

pkg_include_HEADERS = \
$(top_srcdir)/src/navier_stokes/INSProblemCoefs.h \
$(top_srcdir)/src/utilities/ibamr_enums.h \
$(top_srcdir)/src/utilities/ibamr_utilities.h \
$(top_srcdir)/src/utilities/app_namespaces.h \
$(top_srcdir)/src/utilities/namespaces.h

if LIBMESH_ENABLED
pkg_include_HEADERS += \
$(top_srcdir)/src/solid_mechanics/ExplicitFEMechanicsSolver.h
endif

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
$(top_srcdir)/src/IB/GeneralizedIBMethod.C \
$(top_srcdir)/src/IB/IBAnchorPointSpec.C \
$(top_srcdir)/src/IB/IBAnchorPointSpecFactory.C \
$(top_srcdir)/src/IB/IBBeamForceSpec.C \
$(top_srcdir)/src/IB/IBBeamForceSpecFactory.C \
$(top_srcdir)/src/IB/IBEulerianForceFunction.C \
$(top_srcdir)/src/IB/IBEulerianSourceFunction.C \
$(top_srcdir)/src/IB/IBExplicitHierarchyIntegrator.C \
$(top_srcdir)/src/IB/IBHDF5Initializer.C \
$(top_srcdir)/src/IB/IBHierarchyIntegrator.C \
$(top_srcdir)/src/IB/IBInstrumentPanel.C \
$(top_srcdir)/src/IB/IBInstrumentationSpec.C \
$(top_srcdir)/src/IB/IBInstrumentationSpecFactory.C \
$(top_srcdir)/src/IB/IBKirchhoffRodForceGen.C \
$(top_srcdir)/src/IB/IBLagrangianForceStrategy.C \
$(top_srcdir)/src/IB/IBLagrangianForceStrategySet.C \
$(top_srcdir)/src/IB/IBLagrangianSourceStrategy.C \
$(top_srcdir)/src/IB/IBMethod.C \
$(top_srcdir)/src/IB/IBMethodPostProcessStrategy.C \
$(top_srcdir)/src/IB/IBRodForceSpec.C \
$(top_srcdir)/src/IB/IBRodForceSpecFactory.C \
$(top_srcdir)/src/IB/IBSourceSpec.C \
$(top_srcdir)/src/IB/IBSourceSpecFactory.C \
$(top_srcdir)/src/IB/IBSpringForceSpec.C \
$(top_srcdir)/src/IB/IBSpringForceSpecFactory.C \
$(top_srcdir)/src/IB/IBStandardForceGen.C \
$(top_srcdir)/src/IB/IBStandardSourceGen.C \
$(top_srcdir)/src/IB/IBStandardInitializer.C \
$(top_srcdir)/src/IB/IBStrategy.C \
$(top_srcdir)/src/IB/IBStrategySet.C \
$(top_srcdir)/src/IB/IBTargetPointForceSpec.C \
$(top_srcdir)/src/IB/IBTargetPointForceSpecFactory.C \
$(top_srcdir)/src/IB/PenaltyIBMethod.C \
$(top_srcdir)/src/adv_diff/AdvDiffHierarchyIntegrator.C \
$(top_srcdir)/src/adv_diff/AdvDiffHypPatchOps.C \
$(top_srcdir)/src/advect/AdvectHypPatchOps.C \
$(top_srcdir)/src/advect/GodunovAdvector.C \
$(top_srcdir)/src/navier_stokes/ConvectiveOperator.C \
$(top_srcdir)/src/navier_stokes/INSCollocatedCenteredConvectiveOperator.C \
$(top_srcdir)/src/navier_stokes/INSCollocatedHierarchyIntegrator.C \
$(top_srcdir)/src/navier_stokes/INSCollocatedPPMConvectiveOperator.C \
$(top_srcdir)/src/navier_stokes/INSHierarchyIntegrator.C \
$(top_srcdir)/src/navier_stokes/INSIntermediateVelocityBcCoef.C \
$(top_srcdir)/src/navier_stokes/INSPETScMatUtilities.C \
$(top_srcdir)/src/navier_stokes/INSPETScVecUtilities.C \
$(top_srcdir)/src/navier_stokes/INSProjectionBcCoef.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredBlockFactorizationPreconditioner.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredBoxRelaxationFACOperator.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredCenteredConvectiveOperator.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredHierarchyIntegrator.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredPETScLevelSolver.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredPPMConvectiveOperator.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredPhysicalBoundaryHelper.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredPressureBcCoef.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredProjectionPreconditioner.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredStokesOperator.C \
$(top_srcdir)/src/navier_stokes/INSStaggeredVelocityBcCoef.C \
$(top_srcdir)/src/navier_stokes/SpongeLayerForceFunction.C

if LIBMESH_ENABLED
DIM_INDEPENDENT_SOURCES += \
$(top_srcdir)/src/IB/IBFEMethod.C
endif

pkg_include_HEADERS += \
$(top_srcdir)/src/IB/GeneralizedIBMethod.h \
$(top_srcdir)/src/IB/IBAnchorPointSpec.I \
$(top_srcdir)/src/IB/IBAnchorPointSpec.h \
$(top_srcdir)/src/IB/IBBeamForceSpec.I \
$(top_srcdir)/src/IB/IBBeamForceSpec.h \
$(top_srcdir)/src/IB/IBExplicitHierarchyIntegrator.h \
$(top_srcdir)/src/IB/IBHDF5Initializer.h \
$(top_srcdir)/src/IB/IBHierarchyIntegrator.h \
$(top_srcdir)/src/IB/IBInstrumentPanel.h \
$(top_srcdir)/src/IB/IBInstrumentationSpec.I \
$(top_srcdir)/src/IB/IBInstrumentationSpec.h \
$(top_srcdir)/src/IB/IBKirchhoffRodForceGen.h \
$(top_srcdir)/src/IB/IBLagrangianForceStrategy.h \
$(top_srcdir)/src/IB/IBLagrangianForceStrategySet.h \
$(top_srcdir)/src/IB/IBLagrangianSourceStrategy.h \
$(top_srcdir)/src/IB/IBMethod.h \
$(top_srcdir)/src/IB/IBMethodPostProcessStrategy.h \
$(top_srcdir)/src/IB/IBRodForceSpec.I \
$(top_srcdir)/src/IB/IBRodForceSpec.h \
$(top_srcdir)/src/IB/IBSourceSpec.I \
$(top_srcdir)/src/IB/IBSourceSpec.h \
$(top_srcdir)/src/IB/IBSpringForceFunctions.h \
$(top_srcdir)/src/IB/IBSpringForceSpec.I \
$(top_srcdir)/src/IB/IBSpringForceSpec.h \
$(top_srcdir)/src/IB/IBStandardForceGen.h \
$(top_srcdir)/src/IB/IBStandardSourceGen.h \
$(top_srcdir)/src/IB/IBStandardInitializer.h \
$(top_srcdir)/src/IB/IBStrategy.h \
$(top_srcdir)/src/IB/IBStrategySet.h \
$(top_srcdir)/src/IB/IBTargetPointForceSpec.I \
$(top_srcdir)/src/IB/IBTargetPointForceSpec.h \
$(top_srcdir)/src/IB/PenaltyIBMethod.h \
$(top_srcdir)/src/adv_diff/AdvDiffHierarchyIntegrator.h \
$(top_srcdir)/src/adv_diff/AdvDiffHypPatchOps.h \
$(top_srcdir)/src/advect/AdvectHypPatchOps.h \
$(top_srcdir)/src/advect/GodunovAdvector.h \
$(top_srcdir)/src/navier_stokes/ConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSCollocatedCenteredConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSCollocatedHierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSCollocatedPPMConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSHierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSIntermediateVelocityBcCoef.h \
$(top_srcdir)/src/navier_stokes/INSPETScMatUtilities.h \
$(top_srcdir)/src/navier_stokes/INSPETScVecUtilities.h \
$(top_srcdir)/src/navier_stokes/INSProjectionBcCoef.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredBlockFactorizationPreconditioner.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredBoxRelaxationFACOperator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredCenteredConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredHierarchyIntegrator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredPETScLevelSolver.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredPPMConvectiveOperator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredPhysicalBoundaryHelper.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredPressureBcCoef.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredProjectionPreconditioner.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredStokesOperator.h \
$(top_srcdir)/src/navier_stokes/INSStaggeredVelocityBcCoef.h \
$(top_srcdir)/src/navier_stokes/SpongeLayerForceFunction.h

if LIBMESH_ENABLED
pkg_include_HEADERS += \
$(top_srcdir)/src/IB/IBFEMethod.h
endif

libIBAMR2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBAMR2d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/adv_diff/fortran/adv_diff_consdiff2d.f \
$(top_builddir)/src/advect/fortran/advect_centered_derivatives2d.f \
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
$(top_builddir)/src/advect/fortran/advect_centered_derivatives3d.f \
$(top_builddir)/src/advect/fortran/advect_detect3d.f \
$(top_builddir)/src/advect/fortran/advect_diff3d.f \
$(top_builddir)/src/advect/fortran/advect_stable3d.f \
$(top_builddir)/src/advect/fortran/godunov_predictors3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_bdryop3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_divsource3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_stabledt3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_derivatives3d.f \
$(top_builddir)/src/navier_stokes/fortran/navier_stokes_staggered_helpers3d.f
