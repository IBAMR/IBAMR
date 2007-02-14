## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = fortran adv_diff advect lagrangian utilities . tests

lib_LIBRARIES = libIBAMR.a libIBAMR2d.a libIBAMR3d.a

## Dimension-independent library
libIBAMR_a_LIBADD = advect/fortran/godunov_helpers.o lagrangian/fortran/lagrangian_delta.o
libIBAMR_a_SOURCES = # intentionally blank

## Dimension-dependent libraries
libIBAMR2d_a_LIBADD = \
adv_diff/fortran/adv_diff_consdiff2d.o \
adv_diff/libIBAMR_adv_diff2d_a-AdvDiffHierarchyIntegrator.o \
adv_diff/libIBAMR_adv_diff2d_a-AdvDiffHypPatchOps.o \
advect/fortran/godunov_predictors2d.o \
advect/fortran/advect_diff2d.o \
advect/fortran/advect_stable2d.o \
advect/fortran/advect_detect2d.o \
advect/libIBAMR_advect2d_a-GodunovAdvector.o \
advect/libIBAMR_advect2d_a-GodunovHypPatchOps.o \
lagrangian/libIBAMR_lagrangian2d_a-LNodeIndex.o \
lagrangian/libIBAMR_lagrangian2d_a-LEInteractor.o \
lagrangian/libIBAMR_lagrangian2d_a-LNodeIndexData.o \
lagrangian/fortran/lagrangian_interaction2d.o \
lagrangian/libIBAMR_lagrangian2d_a-LNodeIndexSet.o \
lagrangian/libIBAMR_lagrangian2d_a-LNodeIndexDataFactory.o \
lagrangian/libIBAMR_lagrangian2d_a-LDataManager.o \
lagrangian/libIBAMR_lagrangian2d_a-LNodeIndexVariable.o \
lagrangian/libIBAMR_lagrangian2d_a-LagSiloDataWriter.o \
lagrangian/libIBAMR_lagrangian2d_a-LNodeLevelData.o \
lagrangian/libIBAMR_lagrangian2d_a-LNodeInitStrategy.o \
utilities/libIBAMR_utilities2d_a-StashableFactory.o \
utilities/libIBAMR_utilities2d_a-StashableManager.o \
utilities/libIBAMR_utilities2d_a-Stashable.o \
utilities/libIBAMR_utilities2d_a-StashableStream.o \
utilities/libIBAMR_utilities2d_a-get_version.o
libIBAMR2d_a_SOURCES = # intentionally blank

#IB/libIBAMR_IB2d_a-IBEulerianForceSetter.o \
#IB/libIBAMR_IB2d_a-IBEulerianSourceSetter.o \
#IB/libIBAMR_IB2d_a-IBHierarchyIntegrator.o \
#IB/libIBAMR_IB2d_a-IBLagrangianForceStrategy.o \
#IB/libIBAMR_IB2d_a-IBLagrangianSourceStrategy.o \
#IB/libIBAMR_IB2d_a-IBStandardForceGen.o \
#IB/libIBAMR_IB2d_a-IBStandardForceSpecFactory.o \
#IB/libIBAMR_IB2d_a-IBStandardForceSpec.o \
#IB/libIBAMR_IB2d_a-IBStandardInitializer.o \
#
#navier_stokes/fortran/navier_stokes_stabledt2d.o \
#navier_stokes/fortran/navier_stokes_divsource2d.o \
#navier_stokes/libIBAMR_navier_stokes2d_a-INSHierarchyIntegrator.o \
#navier_stokes/libIBAMR_navier_stokes2d_a-HierarchyProjector.o \
#

libIBAMR3d_a_LIBADD = \
adv_diff/libIBAMR_adv_diff3d_a-AdvDiffHypPatchOps.o \
adv_diff/libIBAMR_adv_diff3d_a-AdvDiffHierarchyIntegrator.o \
adv_diff/fortran/adv_diff_consdiff3d.o \
advect/libIBAMR_advect3d_a-GodunovHypPatchOps.o \
advect/fortran/godunov_predictors3d.o \
advect/fortran/advect_detect3d.o \
advect/fortran/advect_stable3d.o \
advect/fortran/advect_diff3d.o \
advect/libIBAMR_advect3d_a-GodunovAdvector.o \
lagrangian/libIBAMR_lagrangian3d_a-LEInteractor.o \
lagrangian/libIBAMR_lagrangian3d_a-LNodeIndexVariable.o \
lagrangian/libIBAMR_lagrangian3d_a-LNodeIndex.o \
lagrangian/fortran/lagrangian_interaction3d.o \
lagrangian/libIBAMR_lagrangian3d_a-LNodeInitStrategy.o \
lagrangian/libIBAMR_lagrangian3d_a-LNodeLevelData.o \
lagrangian/libIBAMR_lagrangian3d_a-LNodeIndexData.o \
lagrangian/libIBAMR_lagrangian3d_a-LNodeIndexSet.o \
lagrangian/libIBAMR_lagrangian3d_a-LNodeIndexDataFactory.o \
lagrangian/libIBAMR_lagrangian3d_a-LagSiloDataWriter.o \
lagrangian/libIBAMR_lagrangian3d_a-LDataManager.o \
utilities/libIBAMR_utilities3d_a-get_version.o \
utilities/libIBAMR_utilities3d_a-Stashable.o \
utilities/libIBAMR_utilities3d_a-StashableFactory.o \
utilities/libIBAMR_utilities3d_a-StashableManager.o \
utilities/libIBAMR_utilities3d_a-StashableStream.o
libIBAMR3d_a_SOURCES = # intentionally blank

#IB/libIBAMR_IB3d_a-IBEulerianForceSetter.o \
#IB/libIBAMR_IB3d_a-IBEulerianSourceSetter.o \
#IB/libIBAMR_IB3d_a-IBHierarchyIntegrator.o \
#IB/libIBAMR_IB3d_a-IBLagrangianForceStrategy.o \
#IB/libIBAMR_IB3d_a-IBLagrangianSourceStrategy.o \
#IB/libIBAMR_IB3d_a-IBStandardForceGen.o \
#IB/libIBAMR_IB3d_a-IBStandardForceSpecFactory.o \
#IB/libIBAMR_IB3d_a-IBStandardForceSpec.o \
#IB/libIBAMR_IB3d_a-IBStandardInitializer.o \
#
#navier_stokes/libIBAMR_navier_stokes3d_a-HierarchyProjector.o \
#navier_stokes/libIBAMR_navier_stokes3d_a-INSHierarchyIntegrator.o \
#navier_stokes/fortran/navier_stokes_stabledt3d.o \
#navier_stokes/fortran/navier_stokes_divsource3d.o \
#
