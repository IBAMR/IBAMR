## Process this file with automake to produce Makefile.in
include $(top_srcdir)/config/Make-rules
SUBDIRS = doc scripts third_party src lib examples
ACLOCAL_AMFLAGS = -I m4

## Standard make targets.
lib: all
examples: lib
	@(cd examples && $(MAKE) $(AM_MAKEFLAGS) $@) || exit 1;
install-exec-local:
	@echo "===================================================================="
	@echo "===================================================================="
	@echo "WARNING: make install is not officially supported by IBTK!"
	@echo "         use at your own risk!"
	@echo "===================================================================="
	@echo "===================================================================="

## Add various "extra" files to the IBTK distribution.
EXTRA_DIST  = $(top_srcdir)/COPYRIGHT
EXTRA_DIST += $(top_srcdir)/autogen.sh
EXTRA_DIST += $(top_srcdir)/doc/Doxyfile.in
EXTRA_DIST += $(top_srcdir)/doc/footer.html
EXTRA_DIST += $(top_srcdir)/doc/header.html
EXTRA_DIST += $(top_srcdir)/doc/libstdc++.tag
EXTRA_DIST += $(top_srcdir)/doc/mainpage.dox
EXTRA_DIST += $(top_srcdir)/doc/online_Doxyfile.in
EXTRA_DIST += $(top_srcdir)/doc/samrai.tag
EXTRA_DIST += $(top_srcdir)/doc/stylesheet.css
EXTRA_DIST += $(top_srcdir)/include/ibtk/BGaussSeidelPreconditioner.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/BGaussSeidelPreconditioner.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/BJacobiPreconditioner.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/BJacobiPreconditioner.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCDivGradHypreLevelSolver.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCDivGradHypreLevelSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCDivGradOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCLaplaceOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCLaplaceOperator2.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCPoissonFACOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCPoissonHypreLevelSolver.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCPoissonHypreLevelSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCPoissonPETScLevelSolver.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/CCPoissonPETScLevelSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartCellDoubleBoundsPreservingConservativeLinearRefine.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartCellDoubleCubicCoarsen.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartCellDoubleQuadraticCFInterpolation.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartCellDoubleQuadraticRefine.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartCellRobinPhysBdryOp.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartExtrapPhysBdryOp.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartGridFunction.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartSideDoubleCubicCoarsen.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartSideDoubleDivPreservingRefine.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartSideDoubleQuadraticCFInterpolation.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartSideDoubleSpecializedLinearRefine.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CartSideRobinPhysBdryOp.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CellNoCornersFillPattern.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CoarseFineBoundaryRefinePatchStrategy.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/CoarsenPatchStrategySet.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/CoarsenPatchStrategySet.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/DebuggingUtilities.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/DiagonalScalingOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/DiagonalScalingPreconditioner.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/DiagonalScalingPreconditioner.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/EdgeDataSynchronization.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/EdgeSynchCopyFillPattern.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/ExtendedRobinBcCoefStrategy.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/FACPreconditioner.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/FACPreconditioner.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/FACPreconditionerStrategy.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/FEDataManager.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/FaceDataSynchronization.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/FaceSynchCopyFillPattern.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/FixedSizedStream.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/FixedSizedStream.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/GeneralOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/HierarchyGhostCellInterpolation.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/HierarchyMathOps.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/HierarchyVecCellDataOpsReal.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/IBTK_CHKERRQ.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/IndexSet.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/IndexSet.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/IndexUtilities.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/IndexUtilities.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/Interpolator.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/Interpolator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/JacobianOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/KrylovLinearSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LDataManager.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/LDataManager.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LEInteractor.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndex.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndex.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexData.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexData.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexDataFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexIterator.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexIterator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexSet.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexSet.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexTransaction.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeIndexVariable.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeInitStrategy.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeLevelData.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/LNodeLevelData.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LagM3DDataWriter.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LagMarker.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/LagMarker.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LagMarkerCoarsen.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LagMarkerRefine.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LagMarkerUtilities.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LagSiloDataWriter.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LinearOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/LinearSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/NewtonKrylovSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/NodeDataSynchronization.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/NodeSynchCopyFillPattern.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/NormOps.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScKrylovLinearSolver.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScKrylovLinearSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScMFFDJacobianOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScMatLOWrapper.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScMatLOWrapper.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScMatOps.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScMatUtilities.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScMultiVec.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScNewtonKrylovSolver.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScNewtonKrylovSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScPCLSWrapper.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScPCLSWrapper.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScSAMRAIVectorReal.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScSAMRAIVectorReal.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScSNESFunctionGOWrapper.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScSNESFunctionGOWrapper.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScSNESJacobianJOWrapper.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScSNESJacobianJOWrapper.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScVecOps.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PETScVecUtilities.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/ParallelEdgeMap.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/ParallelMap.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/ParallelSet.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PatchMathOps.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PatchVecCellDataBasicOps.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PatchVecCellDataOpsReal.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/PhysicalBoundaryUtilities.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/RefinePatchStrategySet.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/RefinePatchStrategySet.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/SCLaplaceOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/SCPoissonFACOperator.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/SCPoissonHypreLevelSolver.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/SCPoissonHypreLevelSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/SCPoissonPETScLevelSolver.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/SCPoissonPETScLevelSolver.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/SideDataSynchronization.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/SideNoCornersFillPattern.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/SideSynchCopyFillPattern.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/StandardTagAndInitStrategySet.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/StandardTagAndInitStrategySet.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/Stashable.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/StashableFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/StashableManager.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/Streamable.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/StreamableFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/StreamableManager.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/StreamableManager.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/TimeDependentLocationIndexRobinBcCoefs.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/VecCellCoarsenAdapter.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/VecCellData.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/VecCellData.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/VecCellDataFactory.I
EXTRA_DIST += $(top_srcdir)/include/ibtk/VecCellDataFactory.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/VecCellRefineAdapter.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/VecCellVariable.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/compiler_hints.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/efficient_add_or_update.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/ibtk_utilities.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/libmesh_utilities.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParser.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserBase.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserBytecode.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserCallback.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserCartGridFunction.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserDef.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserError.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserFixes.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserRobinBcCoefs.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserStack.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserToken.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/muParserTokenReader.h
EXTRA_DIST += $(top_srcdir)/include/ibtk/namespaces.h

## Automatically link header files in $(top_srcdir)/include/ibtk and
## refresh emacs TAGS file.
##
## Note that these targets are for developer use only.
copy_headers:
	cd $(top_srcdir) && scripts/includes.pl
relink_headers:
	cd $(top_srcdir) && scripts/includes.pl --link
refresh_tags:
	cd $(top_srcdir) && find . -name "*.[ChIf]" -print | etags -
