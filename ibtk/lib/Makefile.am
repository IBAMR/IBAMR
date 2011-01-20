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
$(top_srcdir)/src/third_party/muParser/muParser.cpp \
$(top_srcdir)/src/third_party/muParser/muParserBase.cpp \
$(top_srcdir)/src/third_party/muParser/muParserBytecode.cpp \
$(top_srcdir)/src/third_party/muParser/muParserCallback.cpp \
$(top_srcdir)/src/third_party/muParser/muParserError.cpp \
$(top_srcdir)/src/third_party/muParser/muParserTokenReader.cpp \
$(top_srcdir)/src/utilities/Interpolator.C

pkg_include_HEADERS = \
$(top_srcdir)/src/third_party/muParser/muParser.h \
$(top_srcdir)/src/third_party/muParser/muParserBase.h \
$(top_srcdir)/src/third_party/muParser/muParserBytecode.h \
$(top_srcdir)/src/third_party/muParser/muParserCallback.h \
$(top_srcdir)/src/third_party/muParser/muParserDef.h \
$(top_srcdir)/src/third_party/muParser/muParserError.h \
$(top_srcdir)/src/third_party/muParser/muParserFixes.h \
$(top_srcdir)/src/third_party/muParser/muParserStack.h \
$(top_srcdir)/src/third_party/muParser/muParserToken.h \
$(top_srcdir)/src/third_party/muParser/muParserTokenReader.h \
$(top_srcdir)/src/utilities/IBTK_CHKERRQ.h \
$(top_srcdir)/src/utilities/Interpolator.I \
$(top_srcdir)/src/utilities/Interpolator.h \
$(top_srcdir)/src/utilities/compiler_hints.h \
$(top_srcdir)/src/utilities/efficient_add_or_update.h \
$(top_srcdir)/src/utilities/namespaces.h

## Dimension-dependent libraries
DIM_INDEPENDENT_SOURCES = \
$(top_srcdir)/src/boundary/HierarchyGhostCellInterpolation.C \
$(top_srcdir)/src/boundary/cf_interface/CartCellDoubleQuadraticCFInterpolation.C \
$(top_srcdir)/src/boundary/cf_interface/CartSideDoubleQuadraticCFInterpolation.C \
$(top_srcdir)/src/boundary/cf_interface/CoarseFineBoundaryRefinePatchStrategy.C \
$(top_srcdir)/src/boundary/physical_boundary/CartCellRobinPhysBdryOp.C \
$(top_srcdir)/src/boundary/physical_boundary/CartExtrapPhysBdryOp.C \
$(top_srcdir)/src/boundary/physical_boundary/CartSideRobinPhysBdryOp.C \
$(top_srcdir)/src/boundary/physical_boundary/ExtendedRobinBcCoefStrategy.C \
$(top_srcdir)/src/boundary/physical_boundary/PhysicalBoundaryUtilities.C \
$(top_srcdir)/src/boundary/physical_boundary/TimeDependentLocationIndexRobinBcCoefs.C \
$(top_srcdir)/src/boundary/physical_boundary/muParserRobinBcCoefs.C \
$(top_srcdir)/src/coarsen_ops/CartCellDoubleCubicCoarsen.C \
$(top_srcdir)/src/coarsen_ops/CartSideDoubleCubicCoarsen.C \
$(top_srcdir)/src/coarsen_ops/LagMarkerCoarsen.C \
$(top_srcdir)/src/coarsen_ops/VecCellCoarsenAdapter.C \
$(top_srcdir)/src/lagrangian/LDataManager.C \
$(top_srcdir)/src/lagrangian/LEInteractor.C \
$(top_srcdir)/src/lagrangian/LNodeIndex.C \
$(top_srcdir)/src/lagrangian/LNodeIndexData.C \
$(top_srcdir)/src/lagrangian/LNodeIndexDataFactory.C \
$(top_srcdir)/src/lagrangian/LNodeIndexIterator.C \
$(top_srcdir)/src/lagrangian/LNodeIndexSet.C \
$(top_srcdir)/src/lagrangian/LNodeIndexTransaction.C \
$(top_srcdir)/src/lagrangian/LNodeIndexVariable.C \
$(top_srcdir)/src/lagrangian/LNodeInitStrategy.C \
$(top_srcdir)/src/lagrangian/LNodeLevelData.C \
$(top_srcdir)/src/lagrangian/LagM3DDataWriter.C \
$(top_srcdir)/src/lagrangian/LagMarker.C \
$(top_srcdir)/src/lagrangian/LagSiloDataWriter.C \
$(top_srcdir)/src/math/HierarchyMathOps.C \
$(top_srcdir)/src/math/HierarchyVecCellDataOpsReal.C \
$(top_srcdir)/src/math/PETScMatUtilities.C \
$(top_srcdir)/src/math/PETScVecUtilities.C \
$(top_srcdir)/src/math/PatchMathOps.C \
$(top_srcdir)/src/math/PatchVecCellDataBasicOps.C \
$(top_srcdir)/src/math/PatchVecCellDataOpsReal.C \
$(top_srcdir)/src/patch_data/VecCellData.C \
$(top_srcdir)/src/patch_data/VecCellDataFactory.C \
$(top_srcdir)/src/patch_data/VecCellVariable.C \
$(top_srcdir)/src/refine_ops/CartCellDoubleBoundsPreservingConservativeLinearRefine.C \
$(top_srcdir)/src/refine_ops/CartCellDoubleQuadraticRefine.C \
$(top_srcdir)/src/refine_ops/CartSideDoubleDivPreservingRefine.C \
$(top_srcdir)/src/refine_ops/VecCellRefineAdapter.C \
$(top_srcdir)/src/refine_ops/LagMarkerRefine.C \
$(top_srcdir)/src/solvers/impls/BGaussSeidelPreconditioner.C \
$(top_srcdir)/src/solvers/impls/BJacobiPreconditioner.C \
$(top_srcdir)/src/solvers/impls/CCDivGradHypreLevelSolver.C \
$(top_srcdir)/src/solvers/impls/CCDivGradOperator.C \
$(top_srcdir)/src/solvers/impls/CCLaplaceOperator.C \
$(top_srcdir)/src/solvers/impls/CCLaplaceOperator2.C \
$(top_srcdir)/src/solvers/impls/CCPoissonFACOperator.C \
$(top_srcdir)/src/solvers/impls/CCPoissonHypreLevelSolver.C \
$(top_srcdir)/src/solvers/impls/CCPoissonPETScLevelSolver.C \
$(top_srcdir)/src/solvers/impls/DiagonalScalingOperator.C \
$(top_srcdir)/src/solvers/impls/DiagonalScalingPreconditioner.C \
$(top_srcdir)/src/solvers/impls/FACPreconditioner.C \
$(top_srcdir)/src/solvers/impls/PETScKrylovLinearSolver.C \
$(top_srcdir)/src/solvers/impls/PETScMFFDJacobianOperator.C \
$(top_srcdir)/src/solvers/impls/PETScMultiVec.C \
$(top_srcdir)/src/solvers/impls/PETScNewtonKrylovSolver.C \
$(top_srcdir)/src/solvers/impls/SCLaplaceOperator.C \
$(top_srcdir)/src/solvers/impls/SCPoissonFACOperator.C \
$(top_srcdir)/src/solvers/impls/SCPoissonHypreLevelSolver.C \
$(top_srcdir)/src/solvers/impls/SCPoissonPETScLevelSolver.C \
$(top_srcdir)/src/solvers/interfaces/FACPreconditionerStrategy.C \
$(top_srcdir)/src/solvers/interfaces/GeneralOperator.C \
$(top_srcdir)/src/solvers/interfaces/JacobianOperator.C \
$(top_srcdir)/src/solvers/interfaces/KrylovLinearSolver.C \
$(top_srcdir)/src/solvers/interfaces/LinearOperator.C \
$(top_srcdir)/src/solvers/interfaces/LinearSolver.C \
$(top_srcdir)/src/solvers/interfaces/NewtonKrylovSolver.C \
$(top_srcdir)/src/solvers/wrappers/PETScMatLOWrapper.C \
$(top_srcdir)/src/solvers/wrappers/PETScPCLSWrapper.C \
$(top_srcdir)/src/solvers/wrappers/PETScSAMRAIVectorReal.C \
$(top_srcdir)/src/solvers/wrappers/PETScSNESFunctionGOWrapper.C \
$(top_srcdir)/src/solvers/wrappers/PETScSNESJacobianJOWrapper.C \
$(top_srcdir)/src/utilities/CartGridFunction.C \
$(top_srcdir)/src/utilities/CellNoCornersFillPattern.C \
$(top_srcdir)/src/utilities/CoarsenPatchStrategySet.C \
$(top_srcdir)/src/utilities/DebuggingUtilities.C \
$(top_srcdir)/src/utilities/FixedSizedStream.C \
$(top_srcdir)/src/utilities/IndexUtilities.C \
$(top_srcdir)/src/utilities/IndexSet.C \
$(top_srcdir)/src/utilities/LagMarkerUtilities.C \
$(top_srcdir)/src/utilities/NormOps.C \
$(top_srcdir)/src/utilities/PETScMatOps.C \
$(top_srcdir)/src/utilities/PETScVecOps.C \
$(top_srcdir)/src/utilities/ParallelMap.C \
$(top_srcdir)/src/utilities/RefinePatchStrategySet.C \
$(top_srcdir)/src/utilities/SideDataSynchronization.C \
$(top_srcdir)/src/utilities/SideNoCornersFillPattern.C \
$(top_srcdir)/src/utilities/SideSynchCopyFillPattern.C \
$(top_srcdir)/src/utilities/SideSynchCopyTransaction.C \
$(top_srcdir)/src/utilities/SideSynchCopyTransactionFactory.C \
$(top_srcdir)/src/utilities/StandardTagAndInitStrategySet.C \
$(top_srcdir)/src/utilities/Streamable.C \
$(top_srcdir)/src/utilities/StreamableFactory.C \
$(top_srcdir)/src/utilities/StreamableManager.C \
$(top_srcdir)/src/utilities/muParserCartGridFunction.C

if LIBMESH_ENABLED
DIM_INDEPENDENT_SOURCES += \
$(top_srcdir)/src/lagrangian/FEDataManager.C
endif

pkg_include_HEADERS += \
$(top_srcdir)/src/boundary/HierarchyGhostCellInterpolation.h \
$(top_srcdir)/src/boundary/cf_interface/CartCellDoubleQuadraticCFInterpolation.h \
$(top_srcdir)/src/boundary/cf_interface/CartSideDoubleQuadraticCFInterpolation.h \
$(top_srcdir)/src/boundary/cf_interface/CoarseFineBoundaryRefinePatchStrategy.h \
$(top_srcdir)/src/boundary/physical_boundary/CartCellRobinPhysBdryOp.h \
$(top_srcdir)/src/boundary/physical_boundary/CartExtrapPhysBdryOp.h \
$(top_srcdir)/src/boundary/physical_boundary/CartSideRobinPhysBdryOp.h \
$(top_srcdir)/src/boundary/physical_boundary/ExtendedRobinBcCoefStrategy.h \
$(top_srcdir)/src/boundary/physical_boundary/PhysicalBoundaryUtilities.h \
$(top_srcdir)/src/boundary/physical_boundary/TimeDependentLocationIndexRobinBcCoefs.h \
$(top_srcdir)/src/boundary/physical_boundary/muParserRobinBcCoefs.h \
$(top_srcdir)/src/coarsen_ops/CartCellDoubleCubicCoarsen.h \
$(top_srcdir)/src/coarsen_ops/CartSideDoubleCubicCoarsen.h \
$(top_srcdir)/src/coarsen_ops/LagMarkerCoarsen.h \
$(top_srcdir)/src/coarsen_ops/VecCellCoarsenAdapter.h \
$(top_srcdir)/src/lagrangian/LDataManager.h \
$(top_srcdir)/src/lagrangian/LEInteractor.h \
$(top_srcdir)/src/lagrangian/LNodeIndex.I \
$(top_srcdir)/src/lagrangian/LNodeIndex.h \
$(top_srcdir)/src/lagrangian/LNodeIndexData.I \
$(top_srcdir)/src/lagrangian/LNodeIndexData.h \
$(top_srcdir)/src/lagrangian/LNodeIndexDataFactory.h \
$(top_srcdir)/src/lagrangian/LNodeIndexIterator.I \
$(top_srcdir)/src/lagrangian/LNodeIndexIterator.h \
$(top_srcdir)/src/lagrangian/LNodeIndexSet.I \
$(top_srcdir)/src/lagrangian/LNodeIndexSet.h \
$(top_srcdir)/src/lagrangian/LNodeIndexTransaction.h \
$(top_srcdir)/src/lagrangian/LNodeIndexVariable.h \
$(top_srcdir)/src/lagrangian/LNodeInitStrategy.h \
$(top_srcdir)/src/lagrangian/LNodeLevelData.I \
$(top_srcdir)/src/lagrangian/LNodeLevelData.h \
$(top_srcdir)/src/lagrangian/LagM3DDataWriter.h \
$(top_srcdir)/src/lagrangian/LagMarker.I \
$(top_srcdir)/src/lagrangian/LagMarker.h \
$(top_srcdir)/src/lagrangian/LagSiloDataWriter.h \
$(top_srcdir)/src/math/HierarchyMathOps.h \
$(top_srcdir)/src/math/HierarchyVecCellDataOpsReal.h \
$(top_srcdir)/src/math/PETScMatUtilities.h \
$(top_srcdir)/src/math/PETScVecUtilities.h \
$(top_srcdir)/src/math/PatchMathOps.h \
$(top_srcdir)/src/math/PatchVecCellDataBasicOps.h \
$(top_srcdir)/src/math/PatchVecCellDataOpsReal.h \
$(top_srcdir)/src/patch_data/VecCellData.I \
$(top_srcdir)/src/patch_data/VecCellData.h \
$(top_srcdir)/src/patch_data/VecCellDataFactory.I \
$(top_srcdir)/src/patch_data/VecCellDataFactory.h \
$(top_srcdir)/src/patch_data/VecCellVariable.h \
$(top_srcdir)/src/refine_ops/CartCellDoubleBoundsPreservingConservativeLinearRefine.h \
$(top_srcdir)/src/refine_ops/CartCellDoubleQuadraticRefine.h \
$(top_srcdir)/src/refine_ops/CartSideDoubleDivPreservingRefine.h \
$(top_srcdir)/src/refine_ops/VecCellRefineAdapter.h \
$(top_srcdir)/src/refine_ops/LagMarkerRefine.h \
$(top_srcdir)/src/solvers/impls/BGaussSeidelPreconditioner.I \
$(top_srcdir)/src/solvers/impls/BGaussSeidelPreconditioner.h \
$(top_srcdir)/src/solvers/impls/BJacobiPreconditioner.I \
$(top_srcdir)/src/solvers/impls/BJacobiPreconditioner.h \
$(top_srcdir)/src/solvers/impls/CCDivGradHypreLevelSolver.I \
$(top_srcdir)/src/solvers/impls/CCDivGradHypreLevelSolver.h \
$(top_srcdir)/src/solvers/impls/CCDivGradOperator.h \
$(top_srcdir)/src/solvers/impls/CCLaplaceOperator.h \
$(top_srcdir)/src/solvers/impls/CCLaplaceOperator2.h \
$(top_srcdir)/src/solvers/impls/CCPoissonFACOperator.h \
$(top_srcdir)/src/solvers/impls/CCPoissonHypreLevelSolver.I \
$(top_srcdir)/src/solvers/impls/CCPoissonHypreLevelSolver.h \
$(top_srcdir)/src/solvers/impls/CCPoissonPETScLevelSolver.h \
$(top_srcdir)/src/solvers/impls/DiagonalScalingOperator.h \
$(top_srcdir)/src/solvers/impls/DiagonalScalingPreconditioner.I \
$(top_srcdir)/src/solvers/impls/DiagonalScalingPreconditioner.h \
$(top_srcdir)/src/solvers/impls/FACPreconditioner.I \
$(top_srcdir)/src/solvers/impls/FACPreconditioner.h \
$(top_srcdir)/src/solvers/impls/PETScKrylovLinearSolver.I \
$(top_srcdir)/src/solvers/impls/PETScKrylovLinearSolver.h \
$(top_srcdir)/src/solvers/impls/PETScMFFDJacobianOperator.h \
$(top_srcdir)/src/solvers/impls/PETScMultiVec.h \
$(top_srcdir)/src/solvers/impls/PETScNewtonKrylovSolver.I \
$(top_srcdir)/src/solvers/impls/PETScNewtonKrylovSolver.h \
$(top_srcdir)/src/solvers/impls/SCLaplaceOperator.h \
$(top_srcdir)/src/solvers/impls/SCPoissonFACOperator.h \
$(top_srcdir)/src/solvers/impls/SCPoissonHypreLevelSolver.I \
$(top_srcdir)/src/solvers/impls/SCPoissonHypreLevelSolver.h \
$(top_srcdir)/src/solvers/impls/SCPoissonPETScLevelSolver.I \
$(top_srcdir)/src/solvers/impls/SCPoissonPETScLevelSolver.h \
$(top_srcdir)/src/solvers/interfaces/FACPreconditionerStrategy.h \
$(top_srcdir)/src/solvers/interfaces/GeneralOperator.h \
$(top_srcdir)/src/solvers/interfaces/JacobianOperator.h \
$(top_srcdir)/src/solvers/interfaces/KrylovLinearSolver.h \
$(top_srcdir)/src/solvers/interfaces/LinearOperator.h \
$(top_srcdir)/src/solvers/interfaces/LinearSolver.h \
$(top_srcdir)/src/solvers/interfaces/NewtonKrylovSolver.h \
$(top_srcdir)/src/solvers/wrappers/PETScMatLOWrapper.I \
$(top_srcdir)/src/solvers/wrappers/PETScMatLOWrapper.h \
$(top_srcdir)/src/solvers/wrappers/PETScPCLSWrapper.I \
$(top_srcdir)/src/solvers/wrappers/PETScPCLSWrapper.h \
$(top_srcdir)/src/solvers/wrappers/PETScSAMRAIVectorReal.I \
$(top_srcdir)/src/solvers/wrappers/PETScSAMRAIVectorReal.h \
$(top_srcdir)/src/solvers/wrappers/PETScSNESFunctionGOWrapper.I \
$(top_srcdir)/src/solvers/wrappers/PETScSNESFunctionGOWrapper.h \
$(top_srcdir)/src/solvers/wrappers/PETScSNESJacobianJOWrapper.I \
$(top_srcdir)/src/solvers/wrappers/PETScSNESJacobianJOWrapper.h \
$(top_srcdir)/src/utilities/CartGridFunction.h \
$(top_srcdir)/src/utilities/CellNoCornersFillPattern.h \
$(top_srcdir)/src/utilities/CoarsenPatchStrategySet.I \
$(top_srcdir)/src/utilities/CoarsenPatchStrategySet.h \
$(top_srcdir)/src/utilities/DebuggingUtilities.h \
$(top_srcdir)/src/utilities/FixedSizedStream.I \
$(top_srcdir)/src/utilities/FixedSizedStream.h \
$(top_srcdir)/src/utilities/IndexUtilities.I \
$(top_srcdir)/src/utilities/IndexUtilities.h \
$(top_srcdir)/src/utilities/IndexSet.I \
$(top_srcdir)/src/utilities/IndexSet.h \
$(top_srcdir)/src/utilities/LagMarkerUtilities.h \
$(top_srcdir)/src/utilities/NormOps.h \
$(top_srcdir)/src/utilities/PETScMatOps.h \
$(top_srcdir)/src/utilities/PETScVecOps.h \
$(top_srcdir)/src/utilities/ParallelMap.h \
$(top_srcdir)/src/utilities/RefinePatchStrategySet.I \
$(top_srcdir)/src/utilities/RefinePatchStrategySet.h \
$(top_srcdir)/src/utilities/SideDataSynchronization.h \
$(top_srcdir)/src/utilities/SideNoCornersFillPattern.h \
$(top_srcdir)/src/utilities/SideSynchCopyFillPattern.h \
$(top_srcdir)/src/utilities/SideSynchCopyTransaction.h \
$(top_srcdir)/src/utilities/SideSynchCopyTransactionFactory.h \
$(top_srcdir)/src/utilities/StandardTagAndInitStrategySet.I \
$(top_srcdir)/src/utilities/StandardTagAndInitStrategySet.h \
$(top_srcdir)/src/utilities/Streamable.h \
$(top_srcdir)/src/utilities/StreamableFactory.h \
$(top_srcdir)/src/utilities/StreamableManager.I \
$(top_srcdir)/src/utilities/StreamableManager.h \
$(top_srcdir)/src/utilities/muParserCartGridFunction.h

if LIBMESH_ENABLED
DIM_INDEPENDENT_SOURCES += \
$(top_srcdir)/src/lagrangian/FEDataManager.h \
$(top_srcdir)/src/utilities/libmesh_utilities.h
endif

libIBTK2d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=2
libIBTK2d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
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
$(top_builddir)/src/refine_ops/fortran/divpreservingrefine2d.f \
$(top_builddir)/src/solvers/impls/fortran/patchsmoothers2d.f

libIBTK3d_a_CXXFLAGS = $(AM_CXXFLAGS) -DNDIM=3
libIBTK3d_a_SOURCES = $(DIM_INDEPENDENT_SOURCES) \
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
$(top_builddir)/src/refine_ops/fortran/divpreservingrefine3d.f \
$(top_builddir)/src/solvers/impls/fortran/patchsmoothers3d.f
