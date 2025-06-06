## ---------------------------------------------------------------------
##
## Copyright (c) 2011 - 2023 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules

lib_LIBRARIES = libIBTK.a
if SAMRAI2D_ENABLED
lib_LIBRARIES += libIBTK2d.a
endif
if SAMRAI3D_ENABLED
lib_LIBRARIES += libIBTK3d.a
endif

all-local: $(lib_LIBRARIES)
	for file in "lib*.a" ; do \
	  ln -fs $(PWD)/$$file $(PWD)/../src; \
	done

clean-local:
	rm -f $(PWD)/../src/libIBTK.a
	rm -f $(PWD)/../src/libIBTK2d.a
	rm -f $(PWD)/../src/libIBTK3d.a

## Dimension-independent library
libIBTK_a_SOURCES = ../src/dummy.cpp

DIM_INDEPENDENT_SOURCES = $(top_builddir)/src/utilities/fortran/averaging.f

if USING_BUNDLED_MUPARSER
libIBTK_a_SOURCES += \
../contrib/muparser/src/muParserDLL.cpp \
../contrib/muparser/src/muParserTokenReader.cpp \
../contrib/muparser/src/muParserError.cpp \
../contrib/muparser/src/muParserCallback.cpp \
../contrib/muparser/src/muParserBytecode.cpp \
../contrib/muparser/src/muParserBase.cpp \
../contrib/muparser/src/muParserInt.cpp \
../contrib/muparser/src/muParserTest.cpp \
../contrib/muparser/src/muParser.cpp
endif

pkg_include_HEADERS = \
../include/ibtk/IBTK_CHKERRQ.h \
../include/ibtk/app_namespaces.h \
../include/ibtk/compiler_hints.h \
../include/ibtk/ibtk_enums.h \
../include/ibtk/ibtk_utilities.h \
../include/ibtk/namespaces.h

## Dimension-dependent libraries
DIM_DEPENDENT_SOURCES = \
../src/boundary/HierarchyGhostCellInterpolation.cpp \
../src/boundary/cf_interface/CartCellDoubleLinearCFInterpolation.cpp \
../src/boundary/cf_interface/CartCellDoubleQuadraticCFInterpolation.cpp \
../src/boundary/cf_interface/CartSideDoubleQuadraticCFInterpolation.cpp \
../src/boundary/physical_boundary/CartCellRobinPhysBdryOp.cpp \
../src/boundary/physical_boundary/CartExtrapPhysBdryOp.cpp \
../src/boundary/physical_boundary/CartSideRobinPhysBdryOp.cpp \
../src/boundary/physical_boundary/ExtendedRobinBcCoefStrategy.cpp \
../src/boundary/physical_boundary/PhysicalBoundaryUtilities.cpp \
../src/boundary/physical_boundary/RobinPhysBdryPatchStrategy.cpp \
../src/boundary/physical_boundary/StaggeredPhysicalBoundaryHelper.cpp \
../src/boundary/physical_boundary/muParserRobinBcCoefs.cpp \
../src/coarsen_ops/CartCellDoubleCubicCoarsen.cpp \
../src/coarsen_ops/CartSideDoubleCubicCoarsen.cpp \
../src/coarsen_ops/CartSideDoubleRT0Coarsen.cpp \
../src/lagrangian/LData.cpp \
../src/lagrangian/LDataManager.cpp \
../src/lagrangian/LEInteractor.cpp \
../src/lagrangian/LIndexSetData.cpp \
../src/lagrangian/LIndexSetDataFactory.cpp \
../src/lagrangian/LIndexSetVariable.cpp \
../src/lagrangian/LInitStrategy.cpp \
../src/lagrangian/LMesh.cpp \
../src/lagrangian/LNode.cpp \
../src/lagrangian/LNodeIndex.cpp \
../src/lagrangian/LSet.cpp \
../src/lagrangian/LSetData.cpp \
../src/lagrangian/LSetDataFactory.cpp \
../src/lagrangian/LSetDataIterator.cpp \
../src/lagrangian/LSetVariable.cpp \
../src/lagrangian/LSiloDataWriter.cpp \
../src/lagrangian/LTransaction.cpp \
../src/math/HierarchyMathOps.cpp \
../src/math/PETScMatUtilities.cpp \
../src/math/PETScVecUtilities.cpp \
../src/math/PatchMathOps.cpp \
../src/math/PoissonUtilities.cpp \
../src/math/SAMRAIGhostDataAccumulator.cpp \
../src/refine_ops/CartCellDoubleBoundsPreservingConservativeLinearRefine.cpp \
../src/refine_ops/CartCellDoubleQuadraticRefine.cpp \
../src/refine_ops/CartSideDoubleDivPreservingRefine.cpp \
../src/refine_ops/CartSideDoubleRT0Refine.cpp \
../src/refine_ops/CartSideDoubleSpecializedLinearRefine.cpp \
../src/solvers/solver_utilities.cpp \
../src/solvers/impls/BGaussSeidelPreconditioner.cpp \
../src/solvers/impls/BJacobiPreconditioner.cpp \
../src/solvers/impls/CCLaplaceOperator.cpp \
../src/solvers/impls/CCPoissonBoxRelaxationFACOperator.cpp \
../src/solvers/impls/CCPoissonHypreLevelSolver.cpp \
../src/solvers/impls/CCPoissonLevelRelaxationFACOperator.cpp \
../src/solvers/impls/CCPoissonPETScLevelSolver.cpp \
../src/solvers/impls/CCPoissonPointRelaxationFACOperator.cpp \
../src/solvers/impls/CCPoissonSolverManager.cpp \
../src/solvers/impls/FACPreconditioner.cpp \
../src/solvers/impls/KrylovLinearSolverManager.cpp \
../src/solvers/impls/KrylovLinearSolverPoissonSolverInterface.cpp \
../src/solvers/impls/LaplaceOperator.cpp \
../src/solvers/impls/NewtonKrylovSolverManager.cpp \
../src/solvers/impls/PETScKrylovLinearSolver.cpp \
../src/solvers/impls/PETScKrylovPoissonSolver.cpp \
../src/solvers/impls/PETScLevelSolver.cpp \
../src/solvers/impls/PETScMFFDJacobianOperator.cpp \
../src/solvers/impls/PETScNewtonKrylovSolver.cpp \
../src/solvers/impls/PoissonFACPreconditioner.cpp \
../src/solvers/impls/PoissonFACPreconditionerStrategy.cpp \
../src/solvers/impls/PoissonSolver.cpp \
../src/solvers/impls/SCLaplaceOperator.cpp \
../src/solvers/impls/SCPoissonHypreLevelSolver.cpp \
../src/solvers/impls/SCPoissonPETScLevelSolver.cpp \
../src/solvers/impls/SCPoissonPointRelaxationFACOperator.cpp \
../src/solvers/impls/SCPoissonSolverManager.cpp \
../src/solvers/impls/VCSCViscousOpPointRelaxationFACOperator.cpp \
../src/solvers/impls/VCSCViscousOperator.cpp \
../src/solvers/impls/VCSCViscousPETScLevelSolver.cpp \
../src/solvers/interfaces/FACPreconditionerStrategy.cpp \
../src/solvers/interfaces/GeneralOperator.cpp \
../src/solvers/interfaces/GeneralSolver.cpp \
../src/solvers/interfaces/JacobianOperator.cpp \
../src/solvers/interfaces/KrylovLinearSolver.cpp \
../src/solvers/interfaces/LinearOperator.cpp \
../src/solvers/interfaces/LinearSolver.cpp \
../src/solvers/interfaces/NewtonKrylovSolver.cpp \
../src/solvers/wrappers/PETScMatLOWrapper.cpp \
../src/solvers/wrappers/PETScPCLSWrapper.cpp \
../src/solvers/wrappers/PETScSAMRAIVectorReal.cpp \
../src/solvers/wrappers/PETScSNESFunctionGOWrapper.cpp \
../src/solvers/wrappers/PETScSNESJacobianJOWrapper.cpp \
../src/utilities/AppInitializer.cpp \
../src/utilities/CartGridFunction.cpp \
../src/utilities/CartGridFunctionSet.cpp \
../src/utilities/CellNoCornersFillPattern.cpp \
../src/utilities/CoarsenPatchStrategySet.cpp \
../src/utilities/CopyToRootSchedule.cpp \
../src/utilities/CopyToRootTransaction.cpp \
../src/utilities/DebuggingUtilities.cpp \
../src/utilities/EdgeDataSynchronization.cpp \
../src/utilities/EdgeSynchCopyFillPattern.cpp \
../src/utilities/FaceDataSynchronization.cpp \
../src/utilities/FaceSynchCopyFillPattern.cpp \
../src/utilities/FixedSizedStream.cpp \
../src/utilities/HierarchyAveragedDataManager.cpp \
../src/utilities/HierarchyIntegrator.cpp \
../src/utilities/IBTK_MPI.cpp \
../src/utilities/IBTKInit.cpp \
../src/utilities/IndexUtilities.cpp \
../src/utilities/MarkerPatchHierarchy.cpp \
../src/utilities/MergingLoadBalancer.cpp \
../src/utilities/NodeDataSynchronization.cpp \
../src/utilities/NodeSynchCopyFillPattern.cpp \
../src/utilities/NormOps.cpp \
../src/utilities/ParallelEdgeMap.cpp \
../src/utilities/ParallelMap.cpp \
../src/utilities/ParallelSet.cpp \
../src/utilities/PartitioningBox.cpp \
../src/utilities/RefinePatchStrategySet.cpp \
../src/utilities/SAMRAIDataCache.cpp \
../src/utilities/SecondaryHierarchy.cpp \
../src/utilities/SideDataSynchronization.cpp \
../src/utilities/SideNoCornersFillPattern.cpp \
../src/utilities/SideSynchCopyFillPattern.cpp \
../src/utilities/SnapshotCache.cpp \
../src/utilities/snapshot_utilities.cpp \
../src/utilities/StandardTagAndInitStrategySet.cpp \
../src/utilities/Streamable.cpp \
../src/utilities/StreamableManager.cpp \
../src/utilities/box_utilities.cpp \
../src/utilities/ibtk_utilities.cpp \
../src/utilities/muParserCartGridFunction.cpp

if LIBMESH_ENABLED
DIM_DEPENDENT_SOURCES += \
../src/lagrangian/BoxPartitioner.cpp \
../src/lagrangian/StableCentroidPartitioner.cpp \
../src/lagrangian/FEDataInterpolation.cpp \
../src/lagrangian/FEDataManager.cpp \
../src/lagrangian/FEMapping.cpp \
../src/lagrangian/FEProjector.cpp \
../src/lagrangian/FEValues.cpp \
../src/lagrangian/FischerGuess.cpp \
../src/utilities/LibMeshSystemIBVectors.cpp \
../src/utilities/LibMeshSystemVectors.cpp \
../src/utilities/libmesh_utilities.cpp
endif

pkg_include_HEADERS += \
../include/ibtk/AppInitializer.h \
../include/ibtk/BGaussSeidelPreconditioner.h \
../include/ibtk/BJacobiPreconditioner.h \
../include/ibtk/CCLaplaceOperator.h \
../include/ibtk/CCPoissonBoxRelaxationFACOperator.h \
../include/ibtk/CCPoissonHypreLevelSolver.h \
../include/ibtk/CCPoissonLevelRelaxationFACOperator.h \
../include/ibtk/CCPoissonPETScLevelSolver.h \
../include/ibtk/CCPoissonPointRelaxationFACOperator.h \
../include/ibtk/CCPoissonSolverManager.h \
../include/ibtk/CartCellDoubleBoundsPreservingConservativeLinearRefine.h \
../include/ibtk/CartCellDoubleCubicCoarsen.h \
../include/ibtk/CartCellDoubleLinearCFInterpolation.h \
../include/ibtk/CartCellDoubleQuadraticCFInterpolation.h \
../include/ibtk/CartCellDoubleQuadraticRefine.h \
../include/ibtk/CartCellRobinPhysBdryOp.h \
../include/ibtk/CartExtrapPhysBdryOp.h \
../include/ibtk/CartGridFunction.h \
../include/ibtk/CartGridFunctionSet.h \
../include/ibtk/CartSideDoubleCubicCoarsen.h \
../include/ibtk/CartSideDoubleDivPreservingRefine.h \
../include/ibtk/CartSideDoubleQuadraticCFInterpolation.h \
../include/ibtk/CartSideDoubleRT0Coarsen.h \
../include/ibtk/CartSideDoubleRT0Refine.h \
../include/ibtk/CartSideDoubleSpecializedLinearRefine.h \
../include/ibtk/CartSideRobinPhysBdryOp.h \
../include/ibtk/CellNoCornersFillPattern.h \
../include/ibtk/CoarseFineBoundaryRefinePatchStrategy.h \
../include/ibtk/CoarsenPatchStrategySet.h \
../include/ibtk/CopyToRootSchedule.h \
../include/ibtk/CopyToRootTransaction.h \
../include/ibtk/DebuggingUtilities.h \
../include/ibtk/EdgeDataSynchronization.h \
../include/ibtk/EdgeSynchCopyFillPattern.h \
../include/ibtk/ExtendedRobinBcCoefStrategy.h \
../include/ibtk/FACPreconditioner.h \
../include/ibtk/FACPreconditionerStrategy.h \
../include/ibtk/FaceDataSynchronization.h \
../include/ibtk/FaceSynchCopyFillPattern.h \
../include/ibtk/FixedSizedStream.h \
../include/ibtk/GeneralOperator.h \
../include/ibtk/GeneralSolver.h \
../include/ibtk/HierarchyGhostCellInterpolation.h \
../include/ibtk/HierarchyIntegrator.h \
../include/ibtk/HierarchyMathOps.h \
../include/ibtk/IBTK_MPI.h \
../include/ibtk/IBTKInit.h \
../include/ibtk/IndexUtilities.h \
../include/ibtk/JacobianOperator.h \
../include/ibtk/KrylovLinearSolver.h \
../include/ibtk/KrylovLinearSolverManager.h \
../include/ibtk/KrylovLinearSolverPoissonSolverInterface.h \
../include/ibtk/LData.h \
../include/ibtk/LDataManager.h \
../include/ibtk/LEInteractor.h \
../include/ibtk/LIndexSetData.h \
../include/ibtk/LIndexSetDataFactory.h \
../include/ibtk/LIndexSetVariable.h \
../include/ibtk/LInitStrategy.h \
../include/ibtk/LMesh.h \
../include/ibtk/LNode.h \
../include/ibtk/LNodeIndex.h \
../include/ibtk/LNodeIndexSet.h \
../include/ibtk/LNodeIndexSetData.h \
../include/ibtk/LNodeIndexSetDataFactory.h \
../include/ibtk/LNodeIndexSetDataIterator.h \
../include/ibtk/LNodeIndexSetVariable.h \
../include/ibtk/LNodeIndexTransaction.h \
../include/ibtk/LNodeSet.h \
../include/ibtk/LNodeSetData.h \
../include/ibtk/LNodeSetDataFactory.h \
../include/ibtk/LNodeSetDataIterator.h \
../include/ibtk/LNodeSetVariable.h \
../include/ibtk/LNodeTransaction.h \
../include/ibtk/LSet.h \
../include/ibtk/LSetData.h \
../include/ibtk/LSetDataFactory.h \
../include/ibtk/LSetDataIterator.h \
../include/ibtk/LSetVariable.h \
../include/ibtk/LSiloDataWriter.h \
../include/ibtk/LTransaction.h \
../include/ibtk/LaplaceOperator.h \
../include/ibtk/LinearOperator.h \
../include/ibtk/LinearSolver.h \
../include/ibtk/MarkerPatchHierarchy.h \
../include/ibtk/NewtonKrylovSolver.h \
../include/ibtk/NewtonKrylovSolverManager.h \
../include/ibtk/NodeDataSynchronization.h \
../include/ibtk/NodeSynchCopyFillPattern.h \
../include/ibtk/NormOps.h \
../include/ibtk/PETScKrylovLinearSolver.h \
../include/ibtk/PETScKrylovPoissonSolver.h \
../include/ibtk/PETScLevelSolver.h \
../include/ibtk/PETScMFFDJacobianOperator.h \
../include/ibtk/PETScMatLOWrapper.h \
../include/ibtk/PETScMatUtilities.h \
../include/ibtk/PETScNewtonKrylovSolver.h \
../include/ibtk/PETScPCLSWrapper.h \
../include/ibtk/PETScSAMRAIVectorReal.h \
../include/ibtk/PETScSNESFunctionGOWrapper.h \
../include/ibtk/PETScSNESJacobianJOWrapper.h \
../include/ibtk/PETScVecUtilities.h \
../include/ibtk/ParallelEdgeMap.h \
../include/ibtk/ParallelMap.h \
../include/ibtk/ParallelSet.h \
../include/ibtk/PartitioningBox.h \
../include/ibtk/PatchMathOps.h \
../include/ibtk/PhysicalBoundaryUtilities.h \
../include/ibtk/PoissonFACPreconditioner.h \
../include/ibtk/PoissonFACPreconditionerStrategy.h \
../include/ibtk/PoissonSolver.h \
../include/ibtk/PoissonUtilities.h \
../include/ibtk/SAMRAIGhostDataAccumulator.h \
../include/ibtk/RefinePatchStrategySet.h \
../include/ibtk/RobinPhysBdryPatchStrategy.h \
../include/ibtk/SAMRAIDataCache.h \
../include/ibtk/SCLaplaceOperator.h \
../include/ibtk/SCPoissonHypreLevelSolver.h \
../include/ibtk/SCPoissonPETScLevelSolver.h \
../include/ibtk/SCPoissonPointRelaxationFACOperator.h \
../include/ibtk/SCPoissonSolverManager.h \
../include/ibtk/SideDataSynchronization.h \
../include/ibtk/SideNoCornersFillPattern.h \
../include/ibtk/SideSynchCopyFillPattern.h \
../include/ibtk/SnapshotCache.h \
../include/ibtk/snapshot_utilities.h \
../include/ibtk/StaggeredPhysicalBoundaryHelper.h \
../include/ibtk/StandardTagAndInitStrategySet.h \
../include/ibtk/Streamable.h \
../include/ibtk/StreamableFactory.h \
../include/ibtk/StreamableManager.h \
../include/ibtk/VCSCViscousOpPointRelaxationFACOperator.h \
../include/ibtk/VCSCViscousOperator.h \
../include/ibtk/VCSCViscousPETScLevelSolver.h \
../include/ibtk/box_utilities.h \
../include/ibtk/muParserCartGridFunction.h \
../include/ibtk/muParserRobinBcCoefs.h \
../include/ibtk/private/FixedSizedStream-inl.h \
../include/ibtk/private/IndexUtilities-inl.h \
../include/ibtk/private/LData-inl.h \
../include/ibtk/private/LDataManager-inl.h \
../include/ibtk/private/LIndexSetData-inl.h \
../include/ibtk/private/LMesh-inl.h \
../include/ibtk/private/LNode-inl.h \
../include/ibtk/private/LNodeIndex-inl.h \
../include/ibtk/private/LSet-inl.h \
../include/ibtk/private/LSetData-inl.h \
../include/ibtk/private/LSetDataIterator-inl.h \
../include/ibtk/private/PETScSAMRAIVectorReal-inl.h \
../include/ibtk/private/StreamableManager-inl.h

if LIBMESH_ENABLED
DIM_DEPENDENT_SOURCES += \
../include/lagrangian/BoxPartitioner.h \
../include/lagrangian/StableCentroidPartitioner.h \
../include/lagrangian/FEMapping.h \
../include/lagrangian/FEMappingCache.h \
../include/ibtk/FEDataInterpolation.h \
../include/ibtk/FEDataManager.h \
../include/ibtk/FEProjector.h \
../include/ibtk/FEValues.h \
../include/ibtk/LibMeshSystemIBVectors.h \
../include/ibtk/LibMeshSystemVectors.h \
../include/ibtk/libmesh_utilities.h
endif

libIBTK2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBTK2d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) $(DIM_DEPENDENT_SOURCES) \
$(top_builddir)/src/boundary/cf_interface/fortran/linearcfinterpolation2d.f \
$(top_builddir)/src/boundary/cf_interface/fortran/quadcfinterpolation2d.f \
$(top_builddir)/src/boundary/physical_boundary/fortran/cartphysbdryop2d.f \
$(top_builddir)/src/coarsen_ops/fortran/cubiccoarsen2d.f \
$(top_builddir)/src/coarsen_ops/fortran/rt0coarsen2d.f \
$(top_builddir)/src/lagrangian/fortran/lagrangian_interaction2d.f \
$(top_builddir)/src/math/fortran/curl2d.f \
$(top_builddir)/src/math/fortran/div2d.f \
$(top_builddir)/src/math/fortran/flux2d.f \
$(top_builddir)/src/math/fortran/grad2d.f \
$(top_builddir)/src/math/fortran/graddetect2d.f \
$(top_builddir)/src/math/fortran/interp2d.f \
$(top_builddir)/src/math/fortran/laplace2d.f \
$(top_builddir)/src/math/fortran/miscmath2d.f \
$(top_builddir)/src/math/fortran/vclaplace2d.f \
$(top_builddir)/src/math/fortran/rot2d.f \
$(top_builddir)/src/math/fortran/strain2d.f \
$(top_builddir)/src/refine_ops/fortran/cart_side_refine2d.f \
$(top_builddir)/src/refine_ops/fortran/divpreservingrefine2d.f \
$(top_builddir)/src/solvers/impls/fortran/patchsmoothers2d.f

libIBTK3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBTK3d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) $(DIM_DEPENDENT_SOURCES) \
$(top_builddir)/src/boundary/cf_interface/fortran/linearcfinterpolation3d.f \
$(top_builddir)/src/boundary/cf_interface/fortran/quadcfinterpolation3d.f \
$(top_builddir)/src/boundary/physical_boundary/fortran/cartphysbdryop3d.f \
$(top_builddir)/src/coarsen_ops/fortran/cubiccoarsen3d.f \
$(top_builddir)/src/coarsen_ops/fortran/rt0coarsen3d.f \
$(top_builddir)/src/lagrangian/fortran/lagrangian_interaction3d.f \
$(top_builddir)/src/math/fortran/curl3d.f \
$(top_builddir)/src/math/fortran/div3d.f \
$(top_builddir)/src/math/fortran/flux3d.f \
$(top_builddir)/src/math/fortran/grad3d.f \
$(top_builddir)/src/math/fortran/graddetect3d.f \
$(top_builddir)/src/math/fortran/interp3d.f \
$(top_builddir)/src/math/fortran/laplace3d.f \
$(top_builddir)/src/math/fortran/rot3d.f \
$(top_builddir)/src/math/fortran/strain3d.f \
$(top_builddir)/src/math/fortran/miscmath3d.f \
$(top_builddir)/src/math/fortran/vclaplace3d.f \
$(top_builddir)/src/refine_ops/fortran/cart_side_refine3d.f \
$(top_builddir)/src/refine_ops/fortran/divpreservingrefine3d.f \
$(top_builddir)/src/solvers/impls/fortran/patchsmoothers3d.f
