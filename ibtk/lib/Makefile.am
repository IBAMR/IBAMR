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
libIBTK_a_SOURCES = \
$(top_builddir)/src/lagrangian/fortran/lagrangian_delta.f \
$(top_builddir)/src/fortran/minmod.f

pkg_include_HEADERS = \
../src/utilities/IBTK_CHKERRQ.h \
../src/utilities/compiler_hints.h \
../src/utilities/efficient_add_or_update.h \
../src/utilities/ibtk_enums.h \
../src/utilities/ibtk_utilities.h \
../src/utilities/app_namespaces.h \
../src/utilities/namespaces.h

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
../src/boundary/HierarchyGhostCellInterpolation.cpp \
../src/boundary/cf_interface/CartCellDoubleLinearCFInterpolation.cpp \
../src/boundary/cf_interface/CartCellDoubleQuadraticCFInterpolation.cpp \
../src/boundary/cf_interface/CartSideDoubleQuadraticCFInterpolation.cpp \
../src/boundary/cf_interface/CoarseFineBoundaryRefinePatchStrategy.cpp \
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
../src/coarsen_ops/LMarkerCoarsen.cpp \
../src/lagrangian/LData.cpp \
../src/lagrangian/LDataManager.cpp \
../src/lagrangian/LEInteractor.cpp \
../src/lagrangian/LIndexSetData.cpp \
../src/lagrangian/LIndexSetDataFactory.cpp \
../src/lagrangian/LIndexSetVariable.cpp \
../src/lagrangian/LInitStrategy.cpp \
../src/lagrangian/LMarker.cpp \
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
../src/refine_ops/CartCellDoubleBoundsPreservingConservativeLinearRefine.cpp \
../src/refine_ops/CartCellDoubleQuadraticRefine.cpp \
../src/refine_ops/CartSideDoubleDivPreservingRefine.cpp \
../src/refine_ops/CartSideDoubleSpecializedConstantRefine.cpp \
../src/refine_ops/CartSideDoubleSpecializedLinearRefine.cpp \
../src/refine_ops/LMarkerRefine.cpp \
../src/solvers/impls/BGaussSeidelPreconditioner.cpp \
../src/solvers/impls/BJacobiPreconditioner.cpp \
../src/solvers/impls/CCLaplaceOperator.cpp \
../src/solvers/impls/CCPoissonHypreLevelSolver.cpp \
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
../src/solvers/impls/PETScMultiVec.cpp \
../src/solvers/impls/PETScNewtonKrylovSolver.cpp \
../src/solvers/impls/PoissonFACPreconditioner.cpp \
../src/solvers/impls/PoissonFACPreconditionerStrategy.cpp \
../src/solvers/impls/PoissonSolver.cpp \
../src/solvers/impls/SCLaplaceOperator.cpp \
../src/solvers/impls/SCPoissonHypreLevelSolver.cpp \
../src/solvers/impls/SCPoissonPETScLevelSolver.cpp \
../src/solvers/impls/SCPoissonPointRelaxationFACOperator.cpp \
../src/solvers/impls/SCPoissonSolverManager.cpp \
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
../src/utilities/HierarchyIntegrator.cpp \
../src/utilities/IndexUtilities.cpp \
../src/utilities/LMarkerUtilities.cpp \
../src/utilities/NodeDataSynchronization.cpp \
../src/utilities/NodeSynchCopyFillPattern.cpp \
../src/utilities/NormOps.cpp \
../src/utilities/ParallelEdgeMap.cpp \
../src/utilities/ParallelMap.cpp \
../src/utilities/ParallelSet.cpp \
../src/utilities/RefinePatchStrategySet.cpp \
../src/utilities/SideDataSynchronization.cpp \
../src/utilities/SideNoCornersFillPattern.cpp \
../src/utilities/SideSynchCopyFillPattern.cpp \
../src/utilities/StandardTagAndInitStrategySet.cpp \
../src/utilities/Streamable.cpp \
../src/utilities/StreamableFactory.cpp \
../src/utilities/StreamableManager.cpp \
../src/utilities/muParserCartGridFunction.cpp

if LIBMESH_ENABLED
DIM_INDEPENDENT_SOURCES += \
../src/lagrangian/FEDataManager.cpp
endif

pkg_include_HEADERS += \
../src/boundary/HierarchyGhostCellInterpolation.h \
../src/boundary/cf_interface/CartCellDoubleLinearCFInterpolation.h \
../src/boundary/cf_interface/CartCellDoubleQuadraticCFInterpolation.h \
../src/boundary/cf_interface/CartSideDoubleQuadraticCFInterpolation.h \
../src/boundary/cf_interface/CoarseFineBoundaryRefinePatchStrategy.h \
../src/boundary/physical_boundary/CartCellRobinPhysBdryOp.h \
../src/boundary/physical_boundary/CartExtrapPhysBdryOp.h \
../src/boundary/physical_boundary/CartSideRobinPhysBdryOp.h \
../src/boundary/physical_boundary/ExtendedRobinBcCoefStrategy.h \
../src/boundary/physical_boundary/PhysicalBoundaryUtilities.h \
../src/boundary/physical_boundary/RobinPhysBdryPatchStrategy.h \
../src/boundary/physical_boundary/StaggeredPhysicalBoundaryHelper.h \
../src/boundary/physical_boundary/muParserRobinBcCoefs.h \
../src/coarsen_ops/CartCellDoubleCubicCoarsen.h \
../src/coarsen_ops/CartSideDoubleCubicCoarsen.h \
../src/coarsen_ops/LMarkerCoarsen.h \
../src/lagrangian/LData-inl.h \
../src/lagrangian/LData.h \
../src/lagrangian/LDataManager-inl.h \
../src/lagrangian/LDataManager.h \
../src/lagrangian/LEInteractor.h \
../src/lagrangian/LIndexSetData-inl.h \
../src/lagrangian/LIndexSetData.h \
../src/lagrangian/LIndexSetDataFactory.h \
../src/lagrangian/LIndexSetVariable.h \
../src/lagrangian/LInitStrategy.h \
../src/lagrangian/LMarker-inl.h \
../src/lagrangian/LMarker.h \
../src/lagrangian/LMarkerSet.h \
../src/lagrangian/LMarkerSetData.h \
../src/lagrangian/LMarkerSetDataFactory.h \
../src/lagrangian/LMarkerSetDataIterator.h \
../src/lagrangian/LMarkerSetVariable.h \
../src/lagrangian/LMesh-inl.h \
../src/lagrangian/LMesh.h \
../src/lagrangian/LNode-inl.h \
../src/lagrangian/LNode.h \
../src/lagrangian/LNodeIndex-inl.h \
../src/lagrangian/LNodeIndex.h \
../src/lagrangian/LNodeIndexSet.h \
../src/lagrangian/LNodeIndexSetData.h \
../src/lagrangian/LNodeIndexSetDataFactory.h \
../src/lagrangian/LNodeIndexSetDataIterator.h \
../src/lagrangian/LNodeIndexSetVariable.h \
../src/lagrangian/LNodeSet.h \
../src/lagrangian/LNodeSetData.h \
../src/lagrangian/LNodeSetDataFactory.h \
../src/lagrangian/LNodeSetDataIterator.h \
../src/lagrangian/LNodeSetVariable.h \
../src/lagrangian/LSet-inl.h \
../src/lagrangian/LSet.h \
../src/lagrangian/LSetData-inl.h \
../src/lagrangian/LSetData.h \
../src/lagrangian/LSetDataFactory.h \
../src/lagrangian/LSetDataIterator-inl.h \
../src/lagrangian/LSetDataIterator.h \
../src/lagrangian/LSetVariable.h \
../src/lagrangian/LSiloDataWriter.h \
../src/lagrangian/LTransaction.h \
../src/math/HierarchyMathOps.h \
../src/math/PETScMatUtilities.h \
../src/math/PETScVecUtilities.h \
../src/math/PatchMathOps.h \
../src/math/PoissonUtilities.h \
../src/refine_ops/CartCellDoubleBoundsPreservingConservativeLinearRefine.h \
../src/refine_ops/CartCellDoubleQuadraticRefine.h \
../src/refine_ops/CartSideDoubleDivPreservingRefine.h \
../src/refine_ops/CartSideDoubleSpecializedConstantRefine.h \
../src/refine_ops/CartSideDoubleSpecializedLinearRefine.h \
../src/refine_ops/LMarkerRefine.h \
../src/solvers/impls/BGaussSeidelPreconditioner.h \
../src/solvers/impls/BJacobiPreconditioner.h \
../src/solvers/impls/CCLaplaceOperator.h \
../src/solvers/impls/CCPoissonHypreLevelSolver.h \
../src/solvers/impls/CCPoissonPETScLevelSolver.h \
../src/solvers/impls/CCPoissonPointRelaxationFACOperator.h \
../src/solvers/impls/CCPoissonSolverManager.h \
../src/solvers/impls/FACPreconditioner.h \
../src/solvers/impls/KrylovLinearSolverManager.h \
../src/solvers/impls/KrylovLinearSolverPoissonSolverInterface.h \
../src/solvers/impls/LaplaceOperator.h \
../src/solvers/impls/NewtonKrylovSolverManager.h \
../src/solvers/impls/PETScKrylovLinearSolver.h \
../src/solvers/impls/PETScKrylovPoissonSolver.h \
../src/solvers/impls/PETScLevelSolver.h \
../src/solvers/impls/PETScMFFDJacobianOperator.h \
../src/solvers/impls/PETScMultiVec.h \
../src/solvers/impls/PETScNewtonKrylovSolver.h \
../src/solvers/impls/PoissonFACPreconditioner.h \
../src/solvers/impls/PoissonFACPreconditionerStrategy.h \
../src/solvers/impls/PoissonSolver.h \
../src/solvers/impls/SCLaplaceOperator.h \
../src/solvers/impls/SCPoissonHypreLevelSolver.h \
../src/solvers/impls/SCPoissonPETScLevelSolver.h \
../src/solvers/impls/SCPoissonPointRelaxationFACOperator.h \
../src/solvers/impls/SCPoissonSolverManager.h \
../src/solvers/interfaces/FACPreconditionerStrategy.h \
../src/solvers/interfaces/GeneralOperator.h \
../src/solvers/interfaces/GeneralSolver.h \
../src/solvers/interfaces/JacobianOperator.h \
../src/solvers/interfaces/KrylovLinearSolver.h \
../src/solvers/interfaces/LinearOperator.h \
../src/solvers/interfaces/LinearSolver.h \
../src/solvers/interfaces/NewtonKrylovSolver.h \
../src/solvers/wrappers/PETScMatLOWrapper.h \
../src/solvers/wrappers/PETScPCLSWrapper.h \
../src/solvers/wrappers/PETScSAMRAIVectorReal-inl.h \
../src/solvers/wrappers/PETScSAMRAIVectorReal.h \
../src/solvers/wrappers/PETScSNESFunctionGOWrapper.h \
../src/solvers/wrappers/PETScSNESJacobianJOWrapper.h \
../src/utilities/AppInitializer.h \
../src/utilities/CartGridFunction.h \
../src/utilities/CartGridFunctionSet.h \
../src/utilities/CellNoCornersFillPattern.h \
../src/utilities/CoarsenPatchStrategySet.h \
../src/utilities/CopyToRootSchedule.h \
../src/utilities/CopyToRootTransaction.h \
../src/utilities/DebuggingUtilities.h \
../src/utilities/EdgeDataSynchronization.h \
../src/utilities/EdgeSynchCopyFillPattern.h \
../src/utilities/FaceDataSynchronization.h \
../src/utilities/FaceSynchCopyFillPattern.h \
../src/utilities/FixedSizedStream-inl.h \
../src/utilities/FixedSizedStream.h \
../src/utilities/HierarchyIntegrator.h \
../src/utilities/IndexUtilities-inl.h \
../src/utilities/IndexUtilities.h \
../src/utilities/LMarkerUtilities.h \
../src/utilities/NodeDataSynchronization.h \
../src/utilities/NodeSynchCopyFillPattern.h \
../src/utilities/NormOps.h \
../src/utilities/ParallelEdgeMap.h \
../src/utilities/ParallelMap.h \
../src/utilities/ParallelSet.h \
../src/utilities/RefinePatchStrategySet.h \
../src/utilities/SideDataSynchronization.h \
../src/utilities/SideNoCornersFillPattern.h \
../src/utilities/SideSynchCopyFillPattern.h \
../src/utilities/StandardTagAndInitStrategySet.h \
../src/utilities/Streamable.h \
../src/utilities/StreamableFactory.h \
../src/utilities/StreamableManager-inl.h \
../src/utilities/StreamableManager.h \
../src/utilities/muParserCartGridFunction.h

if LIBMESH_ENABLED
DIM_INDEPENDENT_SOURCES += \
../src/lagrangian/FEDataManager.h \
../src/utilities/libmesh_utilities.h
endif

libIBTK2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBTK2d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/boundary/cf_interface/fortran/linearcfinterpolation2d.f \
$(top_builddir)/src/boundary/cf_interface/fortran/quadcfinterpolation2d.f \
$(top_builddir)/src/boundary/physical_boundary/fortran/cartphysbdryop2d.f \
$(top_builddir)/src/coarsen_ops/fortran/cubiccoarsen2d.f \
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
$(top_builddir)/src/refine_ops/fortran/cart_side_refine2d.f \
$(top_builddir)/src/refine_ops/fortran/divpreservingrefine2d.f \
$(top_builddir)/src/solvers/impls/fortran/patchsmoothers2d.f

libIBTK3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBTK3d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
$(top_builddir)/src/boundary/cf_interface/fortran/linearcfinterpolation3d.f \
$(top_builddir)/src/boundary/cf_interface/fortran/quadcfinterpolation3d.f \
$(top_builddir)/src/boundary/physical_boundary/fortran/cartphysbdryop3d.f \
$(top_builddir)/src/coarsen_ops/fortran/cubiccoarsen3d.f \
$(top_builddir)/src/lagrangian/fortran/lagrangian_interaction3d.f \
$(top_builddir)/src/math/fortran/curl3d.f \
$(top_builddir)/src/math/fortran/div3d.f \
$(top_builddir)/src/math/fortran/flux3d.f \
$(top_builddir)/src/math/fortran/grad3d.f \
$(top_builddir)/src/math/fortran/graddetect3d.f \
$(top_builddir)/src/math/fortran/interp3d.f \
$(top_builddir)/src/math/fortran/laplace3d.f \
$(top_builddir)/src/math/fortran/rot3d.f \
$(top_builddir)/src/math/fortran/miscmath3d.f \
$(top_builddir)/src/refine_ops/fortran/cart_side_refine3d.f \
$(top_builddir)/src/refine_ops/fortran/divpreservingrefine3d.f \
$(top_builddir)/src/solvers/impls/fortran/patchsmoothers3d.f
