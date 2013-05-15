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
$(top_srcdir)/src/utilities/Interpolator.cpp \
$(top_builddir)/src/fortran/minmod.f

pkg_include_HEADERS = \
$(top_srcdir)/src/utilities/IBTK_CHKERRQ.h \
$(top_srcdir)/src/utilities/Interpolator-inl.h \
$(top_srcdir)/src/utilities/Interpolator.h \
$(top_srcdir)/src/utilities/compiler_hints.h \
$(top_srcdir)/src/utilities/efficient_add_or_update.h \
$(top_srcdir)/src/utilities/ibtk_enums.h \
$(top_srcdir)/src/utilities/ibtk_utilities.h \
$(top_srcdir)/src/utilities/app_namespaces.h \
$(top_srcdir)/src/utilities/namespaces.h

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
$(top_srcdir)/src/boundary/HierarchyGhostCellInterpolation.cpp \
$(top_srcdir)/src/boundary/cf_interface/CartCellDoubleLinearCFInterpolation.cpp \
$(top_srcdir)/src/boundary/cf_interface/CartCellDoubleQuadraticCFInterpolation.cpp \
$(top_srcdir)/src/boundary/cf_interface/CartSideDoubleQuadraticCFInterpolation.cpp \
$(top_srcdir)/src/boundary/cf_interface/CoarseFineBoundaryRefinePatchStrategy.cpp \
$(top_srcdir)/src/boundary/physical_boundary/CartCellRobinPhysBdryOp.cpp \
$(top_srcdir)/src/boundary/physical_boundary/CartExtrapPhysBdryOp.cpp \
$(top_srcdir)/src/boundary/physical_boundary/CartSideRobinPhysBdryOp.cpp \
$(top_srcdir)/src/boundary/physical_boundary/ExtendedRobinBcCoefStrategy.cpp \
$(top_srcdir)/src/boundary/physical_boundary/PhysicalBoundaryUtilities.cpp \
$(top_srcdir)/src/boundary/physical_boundary/RobinPhysBdryPatchStrategy.cpp \
$(top_srcdir)/src/boundary/physical_boundary/StaggeredPhysicalBoundaryHelper.cpp \
$(top_srcdir)/src/boundary/physical_boundary/muParserRobinBcCoefs.cpp \
$(top_srcdir)/src/coarsen_ops/CartCellDoubleCubicCoarsen.cpp \
$(top_srcdir)/src/coarsen_ops/CartSideDoubleCubicCoarsen.cpp \
$(top_srcdir)/src/coarsen_ops/LMarkerCoarsen.cpp \
$(top_srcdir)/src/lagrangian/LData.cpp \
$(top_srcdir)/src/lagrangian/LDataManager.cpp \
$(top_srcdir)/src/lagrangian/LEInteractor.cpp \
$(top_srcdir)/src/lagrangian/LIndexSetData.cpp \
$(top_srcdir)/src/lagrangian/LIndexSetDataFactory.cpp \
$(top_srcdir)/src/lagrangian/LIndexSetVariable.cpp \
$(top_srcdir)/src/lagrangian/LInitStrategy.cpp \
$(top_srcdir)/src/lagrangian/LMarker.cpp \
$(top_srcdir)/src/lagrangian/LMesh.cpp \
$(top_srcdir)/src/lagrangian/LNode.cpp \
$(top_srcdir)/src/lagrangian/LNodeIndex.cpp \
$(top_srcdir)/src/lagrangian/LSet.cpp \
$(top_srcdir)/src/lagrangian/LSetData.cpp \
$(top_srcdir)/src/lagrangian/LSetDataFactory.cpp \
$(top_srcdir)/src/lagrangian/LSetDataIterator.cpp \
$(top_srcdir)/src/lagrangian/LSetVariable.cpp \
$(top_srcdir)/src/lagrangian/LSiloDataWriter.cpp \
$(top_srcdir)/src/lagrangian/LTransaction.cpp \
$(top_srcdir)/src/math/HierarchyMathOps.cpp \
$(top_srcdir)/src/math/PETScMatUtilities.cpp \
$(top_srcdir)/src/math/PETScVecUtilities.cpp \
$(top_srcdir)/src/math/PatchMathOps.cpp \
$(top_srcdir)/src/math/PoissonUtilities.cpp \
$(top_srcdir)/src/refine_ops/CartCellDoubleBoundsPreservingConservativeLinearRefine.cpp \
$(top_srcdir)/src/refine_ops/CartCellDoubleQuadraticRefine.cpp \
$(top_srcdir)/src/refine_ops/CartSideDoubleDivPreservingRefine.cpp \
$(top_srcdir)/src/refine_ops/CartSideDoubleSpecializedConstantRefine.cpp \
$(top_srcdir)/src/refine_ops/CartSideDoubleSpecializedLinearRefine.cpp \
$(top_srcdir)/src/refine_ops/LMarkerRefine.cpp \
$(top_srcdir)/src/solvers/impls/BGaussSeidelPreconditioner.cpp \
$(top_srcdir)/src/solvers/impls/BJacobiPreconditioner.cpp \
$(top_srcdir)/src/solvers/impls/CCLaplaceOperator.cpp \
$(top_srcdir)/src/solvers/impls/CCPoissonHypreLevelSolver.cpp \
$(top_srcdir)/src/solvers/impls/CCPoissonPETScLevelSolver.cpp \
$(top_srcdir)/src/solvers/impls/CCPoissonPointRelaxationFACOperator.cpp \
$(top_srcdir)/src/solvers/impls/CCPoissonSolverManager.cpp \
$(top_srcdir)/src/solvers/impls/FACPreconditioner.cpp \
$(top_srcdir)/src/solvers/impls/KrylovLinearSolverManager.cpp \
$(top_srcdir)/src/solvers/impls/KrylovLinearSolverPoissonSolverInterface.cpp \
$(top_srcdir)/src/solvers/impls/LaplaceOperator.cpp \
$(top_srcdir)/src/solvers/impls/NewtonKrylovSolverManager.cpp \
$(top_srcdir)/src/solvers/impls/PETScKrylovLinearSolver.cpp \
$(top_srcdir)/src/solvers/impls/PETScKrylovPoissonSolver.cpp \
$(top_srcdir)/src/solvers/impls/PETScLevelSolver.cpp \
$(top_srcdir)/src/solvers/impls/PETScMFFDJacobianOperator.cpp \
$(top_srcdir)/src/solvers/impls/PETScMultiVec.cpp \
$(top_srcdir)/src/solvers/impls/PETScNewtonKrylovSolver.cpp \
$(top_srcdir)/src/solvers/impls/PoissonFACPreconditioner.cpp \
$(top_srcdir)/src/solvers/impls/PoissonFACPreconditionerStrategy.cpp \
$(top_srcdir)/src/solvers/impls/PoissonSolver.cpp \
$(top_srcdir)/src/solvers/impls/SCLaplaceOperator.cpp \
$(top_srcdir)/src/solvers/impls/SCPoissonHypreLevelSolver.cpp \
$(top_srcdir)/src/solvers/impls/SCPoissonPETScLevelSolver.cpp \
$(top_srcdir)/src/solvers/impls/SCPoissonPointRelaxationFACOperator.cpp \
$(top_srcdir)/src/solvers/impls/SCPoissonSolverManager.cpp \
$(top_srcdir)/src/solvers/interfaces/FACPreconditionerStrategy.cpp \
$(top_srcdir)/src/solvers/interfaces/GeneralOperator.cpp \
$(top_srcdir)/src/solvers/interfaces/GeneralSolver.cpp \
$(top_srcdir)/src/solvers/interfaces/JacobianOperator.cpp \
$(top_srcdir)/src/solvers/interfaces/KrylovLinearSolver.cpp \
$(top_srcdir)/src/solvers/interfaces/LinearOperator.cpp \
$(top_srcdir)/src/solvers/interfaces/LinearSolver.cpp \
$(top_srcdir)/src/solvers/interfaces/NewtonKrylovSolver.cpp \
$(top_srcdir)/src/solvers/wrappers/PETScMatLOWrapper.cpp \
$(top_srcdir)/src/solvers/wrappers/PETScPCLSWrapper.cpp \
$(top_srcdir)/src/solvers/wrappers/PETScSAMRAIVectorReal.cpp \
$(top_srcdir)/src/solvers/wrappers/PETScSNESFunctionGOWrapper.cpp \
$(top_srcdir)/src/solvers/wrappers/PETScSNESJacobianJOWrapper.cpp \
$(top_srcdir)/src/utilities/AppInitializer.cpp \
$(top_srcdir)/src/utilities/CartGridFunction.cpp \
$(top_srcdir)/src/utilities/CartGridFunctionSet.cpp \
$(top_srcdir)/src/utilities/CellNoCornersFillPattern.cpp \
$(top_srcdir)/src/utilities/CoarsenPatchStrategySet.cpp \
$(top_srcdir)/src/utilities/CopyToRootSchedule.cpp \
$(top_srcdir)/src/utilities/CopyToRootTransaction.cpp \
$(top_srcdir)/src/utilities/DebuggingUtilities.cpp \
$(top_srcdir)/src/utilities/EdgeDataSynchronization.cpp \
$(top_srcdir)/src/utilities/EdgeSynchCopyFillPattern.cpp \
$(top_srcdir)/src/utilities/FaceDataSynchronization.cpp \
$(top_srcdir)/src/utilities/FaceSynchCopyFillPattern.cpp \
$(top_srcdir)/src/utilities/FixedSizedStream.cpp \
$(top_srcdir)/src/utilities/HierarchyIntegrator.cpp \
$(top_srcdir)/src/utilities/IndexUtilities.cpp \
$(top_srcdir)/src/utilities/LMarkerUtilities.cpp \
$(top_srcdir)/src/utilities/NodeDataSynchronization.cpp \
$(top_srcdir)/src/utilities/NodeSynchCopyFillPattern.cpp \
$(top_srcdir)/src/utilities/NormOps.cpp \
$(top_srcdir)/src/utilities/ParallelEdgeMap.cpp \
$(top_srcdir)/src/utilities/ParallelMap.cpp \
$(top_srcdir)/src/utilities/ParallelSet.cpp \
$(top_srcdir)/src/utilities/RefinePatchStrategySet.cpp \
$(top_srcdir)/src/utilities/SideDataSynchronization.cpp \
$(top_srcdir)/src/utilities/SideNoCornersFillPattern.cpp \
$(top_srcdir)/src/utilities/SideSynchCopyFillPattern.cpp \
$(top_srcdir)/src/utilities/StandardTagAndInitStrategySet.cpp \
$(top_srcdir)/src/utilities/Streamable.cpp \
$(top_srcdir)/src/utilities/StreamableFactory.cpp \
$(top_srcdir)/src/utilities/StreamableManager.cpp \
$(top_srcdir)/src/utilities/muParserCartGridFunction.cpp

if LIBMESH_ENABLED
DIM_INDEPENDENT_SOURCES += \
$(top_srcdir)/src/lagrangian/FEDataManager.cpp
endif

pkg_include_HEADERS += \
$(top_srcdir)/src/boundary/HierarchyGhostCellInterpolation.h \
$(top_srcdir)/src/boundary/cf_interface/CartCellDoubleLinearCFInterpolation.h \
$(top_srcdir)/src/boundary/cf_interface/CartCellDoubleQuadraticCFInterpolation.h \
$(top_srcdir)/src/boundary/cf_interface/CartSideDoubleQuadraticCFInterpolation.h \
$(top_srcdir)/src/boundary/cf_interface/CoarseFineBoundaryRefinePatchStrategy.h \
$(top_srcdir)/src/boundary/physical_boundary/CartCellRobinPhysBdryOp.h \
$(top_srcdir)/src/boundary/physical_boundary/CartExtrapPhysBdryOp.h \
$(top_srcdir)/src/boundary/physical_boundary/CartSideRobinPhysBdryOp.h \
$(top_srcdir)/src/boundary/physical_boundary/ExtendedRobinBcCoefStrategy.h \
$(top_srcdir)/src/boundary/physical_boundary/PhysicalBoundaryUtilities.h \
$(top_srcdir)/src/boundary/physical_boundary/RobinPhysBdryPatchStrategy.h \
$(top_srcdir)/src/boundary/physical_boundary/StaggeredPhysicalBoundaryHelper.h \
$(top_srcdir)/src/boundary/physical_boundary/muParserRobinBcCoefs.h \
$(top_srcdir)/src/coarsen_ops/CartCellDoubleCubicCoarsen.h \
$(top_srcdir)/src/coarsen_ops/CartSideDoubleCubicCoarsen.h \
$(top_srcdir)/src/coarsen_ops/LMarkerCoarsen.h \
$(top_srcdir)/src/lagrangian/LData-inl.h \
$(top_srcdir)/src/lagrangian/LData.h \
$(top_srcdir)/src/lagrangian/LDataManager-inl.h \
$(top_srcdir)/src/lagrangian/LDataManager.h \
$(top_srcdir)/src/lagrangian/LEInteractor.h \
$(top_srcdir)/src/lagrangian/LIndexSetData-inl.h \
$(top_srcdir)/src/lagrangian/LIndexSetData.h \
$(top_srcdir)/src/lagrangian/LIndexSetDataFactory.h \
$(top_srcdir)/src/lagrangian/LIndexSetVariable.h \
$(top_srcdir)/src/lagrangian/LInitStrategy.h \
$(top_srcdir)/src/lagrangian/LMarker-inl.h \
$(top_srcdir)/src/lagrangian/LMarker.h \
$(top_srcdir)/src/lagrangian/LMarkerSet.h \
$(top_srcdir)/src/lagrangian/LMarkerSetData.h \
$(top_srcdir)/src/lagrangian/LMarkerSetDataFactory.h \
$(top_srcdir)/src/lagrangian/LMarkerSetDataIterator.h \
$(top_srcdir)/src/lagrangian/LMarkerSetVariable.h \
$(top_srcdir)/src/lagrangian/LMesh-inl.h \
$(top_srcdir)/src/lagrangian/LMesh.h \
$(top_srcdir)/src/lagrangian/LNode-inl.h \
$(top_srcdir)/src/lagrangian/LNode.h \
$(top_srcdir)/src/lagrangian/LNodeIndex-inl.h \
$(top_srcdir)/src/lagrangian/LNodeIndex.h \
$(top_srcdir)/src/lagrangian/LNodeIndexSet.h \
$(top_srcdir)/src/lagrangian/LNodeIndexSetData.h \
$(top_srcdir)/src/lagrangian/LNodeIndexSetDataFactory.h \
$(top_srcdir)/src/lagrangian/LNodeIndexSetDataIterator.h \
$(top_srcdir)/src/lagrangian/LNodeIndexSetVariable.h \
$(top_srcdir)/src/lagrangian/LNodeSet.h \
$(top_srcdir)/src/lagrangian/LNodeSetData.h \
$(top_srcdir)/src/lagrangian/LNodeSetDataFactory.h \
$(top_srcdir)/src/lagrangian/LNodeSetDataIterator.h \
$(top_srcdir)/src/lagrangian/LNodeSetVariable.h \
$(top_srcdir)/src/lagrangian/LSet-inl.h \
$(top_srcdir)/src/lagrangian/LSet.h \
$(top_srcdir)/src/lagrangian/LSetData-inl.h \
$(top_srcdir)/src/lagrangian/LSetData.h \
$(top_srcdir)/src/lagrangian/LSetDataFactory.h \
$(top_srcdir)/src/lagrangian/LSetDataIterator-inl.h \
$(top_srcdir)/src/lagrangian/LSetDataIterator.h \
$(top_srcdir)/src/lagrangian/LSetVariable.h \
$(top_srcdir)/src/lagrangian/LSiloDataWriter.h \
$(top_srcdir)/src/lagrangian/LTransaction.h \
$(top_srcdir)/src/math/HierarchyMathOps.h \
$(top_srcdir)/src/math/PETScMatUtilities.h \
$(top_srcdir)/src/math/PETScVecUtilities.h \
$(top_srcdir)/src/math/PatchMathOps.h \
$(top_srcdir)/src/math/PoissonUtilities.h \
$(top_srcdir)/src/refine_ops/CartCellDoubleBoundsPreservingConservativeLinearRefine.h \
$(top_srcdir)/src/refine_ops/CartCellDoubleQuadraticRefine.h \
$(top_srcdir)/src/refine_ops/CartSideDoubleDivPreservingRefine.h \
$(top_srcdir)/src/refine_ops/CartSideDoubleSpecializedConstantRefine.h \
$(top_srcdir)/src/refine_ops/CartSideDoubleSpecializedLinearRefine.h \
$(top_srcdir)/src/refine_ops/LMarkerRefine.h \
$(top_srcdir)/src/solvers/impls/BGaussSeidelPreconditioner.h \
$(top_srcdir)/src/solvers/impls/BJacobiPreconditioner.h \
$(top_srcdir)/src/solvers/impls/CCLaplaceOperator.h \
$(top_srcdir)/src/solvers/impls/CCPoissonHypreLevelSolver.h \
$(top_srcdir)/src/solvers/impls/CCPoissonPETScLevelSolver.h \
$(top_srcdir)/src/solvers/impls/CCPoissonPointRelaxationFACOperator.h \
$(top_srcdir)/src/solvers/impls/CCPoissonSolverManager.h \
$(top_srcdir)/src/solvers/impls/FACPreconditioner.h \
$(top_srcdir)/src/solvers/impls/KrylovLinearSolverManager.h \
$(top_srcdir)/src/solvers/impls/KrylovLinearSolverPoissonSolverInterface.h \
$(top_srcdir)/src/solvers/impls/LaplaceOperator.h \
$(top_srcdir)/src/solvers/impls/NewtonKrylovSolverManager.h \
$(top_srcdir)/src/solvers/impls/PETScKrylovLinearSolver.h \
$(top_srcdir)/src/solvers/impls/PETScKrylovPoissonSolver.h \
$(top_srcdir)/src/solvers/impls/PETScLevelSolver.h \
$(top_srcdir)/src/solvers/impls/PETScMFFDJacobianOperator.h \
$(top_srcdir)/src/solvers/impls/PETScMultiVec.h \
$(top_srcdir)/src/solvers/impls/PETScNewtonKrylovSolver.h \
$(top_srcdir)/src/solvers/impls/PoissonFACPreconditioner.h \
$(top_srcdir)/src/solvers/impls/PoissonFACPreconditionerStrategy.h \
$(top_srcdir)/src/solvers/impls/PoissonSolver.h \
$(top_srcdir)/src/solvers/impls/SCLaplaceOperator.h \
$(top_srcdir)/src/solvers/impls/SCPoissonHypreLevelSolver.h \
$(top_srcdir)/src/solvers/impls/SCPoissonPETScLevelSolver.h \
$(top_srcdir)/src/solvers/impls/SCPoissonPointRelaxationFACOperator.h \
$(top_srcdir)/src/solvers/impls/SCPoissonSolverManager.h \
$(top_srcdir)/src/solvers/interfaces/FACPreconditionerStrategy.h \
$(top_srcdir)/src/solvers/interfaces/GeneralOperator.h \
$(top_srcdir)/src/solvers/interfaces/GeneralSolver.h \
$(top_srcdir)/src/solvers/interfaces/JacobianOperator.h \
$(top_srcdir)/src/solvers/interfaces/KrylovLinearSolver.h \
$(top_srcdir)/src/solvers/interfaces/LinearOperator.h \
$(top_srcdir)/src/solvers/interfaces/LinearSolver.h \
$(top_srcdir)/src/solvers/interfaces/NewtonKrylovSolver.h \
$(top_srcdir)/src/solvers/wrappers/PETScMatLOWrapper.h \
$(top_srcdir)/src/solvers/wrappers/PETScPCLSWrapper.h \
$(top_srcdir)/src/solvers/wrappers/PETScSAMRAIVectorReal-inl.h \
$(top_srcdir)/src/solvers/wrappers/PETScSAMRAIVectorReal.h \
$(top_srcdir)/src/solvers/wrappers/PETScSNESFunctionGOWrapper.h \
$(top_srcdir)/src/solvers/wrappers/PETScSNESJacobianJOWrapper.h \
$(top_srcdir)/src/utilities/AppInitializer.h \
$(top_srcdir)/src/utilities/CartGridFunction.h \
$(top_srcdir)/src/utilities/CartGridFunctionSet.h \
$(top_srcdir)/src/utilities/CellNoCornersFillPattern.h \
$(top_srcdir)/src/utilities/CoarsenPatchStrategySet.h \
$(top_srcdir)/src/utilities/CopyToRootSchedule.h \
$(top_srcdir)/src/utilities/CopyToRootTransaction.h \
$(top_srcdir)/src/utilities/DebuggingUtilities.h \
$(top_srcdir)/src/utilities/EdgeDataSynchronization.h \
$(top_srcdir)/src/utilities/EdgeSynchCopyFillPattern.h \
$(top_srcdir)/src/utilities/FaceDataSynchronization.h \
$(top_srcdir)/src/utilities/FaceSynchCopyFillPattern.h \
$(top_srcdir)/src/utilities/FixedSizedStream-inl.h \
$(top_srcdir)/src/utilities/FixedSizedStream.h \
$(top_srcdir)/src/utilities/HierarchyIntegrator.h \
$(top_srcdir)/src/utilities/IndexUtilities-inl.h \
$(top_srcdir)/src/utilities/IndexUtilities.h \
$(top_srcdir)/src/utilities/LMarkerUtilities.h \
$(top_srcdir)/src/utilities/NodeDataSynchronization.h \
$(top_srcdir)/src/utilities/NodeSynchCopyFillPattern.h \
$(top_srcdir)/src/utilities/NormOps.h \
$(top_srcdir)/src/utilities/ParallelEdgeMap.h \
$(top_srcdir)/src/utilities/ParallelMap.h \
$(top_srcdir)/src/utilities/ParallelSet.h \
$(top_srcdir)/src/utilities/RefinePatchStrategySet.h \
$(top_srcdir)/src/utilities/SideDataSynchronization.h \
$(top_srcdir)/src/utilities/SideNoCornersFillPattern.h \
$(top_srcdir)/src/utilities/SideSynchCopyFillPattern.h \
$(top_srcdir)/src/utilities/StandardTagAndInitStrategySet.h \
$(top_srcdir)/src/utilities/Streamable.h \
$(top_srcdir)/src/utilities/StreamableFactory.h \
$(top_srcdir)/src/utilities/StreamableManager-inl.h \
$(top_srcdir)/src/utilities/StreamableManager.h \
$(top_srcdir)/src/utilities/muParserCartGridFunction.h

if LIBMESH_ENABLED
DIM_INDEPENDENT_SOURCES += \
$(top_srcdir)/src/lagrangian/FEDataManager.h \
$(top_srcdir)/src/utilities/libmesh_utilities.h
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
$(top_builddir)/src/math/fortran/miscmath3d.f \
$(top_builddir)/src/refine_ops/fortran/cart_side_refine3d.f \
$(top_builddir)/src/refine_ops/fortran/divpreservingrefine3d.f \
$(top_builddir)/src/solvers/impls/fortran/patchsmoothers3d.f
