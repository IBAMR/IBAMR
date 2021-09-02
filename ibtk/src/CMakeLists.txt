## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2021 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

ADD_LIBRARY(IBTK2d)
ADD_LIBRARY(IBTK3d)

# 'base' names (without, 2d.f.m4 or 3d.f.m4) for all the fortran that is
# explicitly compiled
SET(FORTRAN_SRC_BASE
  boundary/cf_interface/fortran/linearcfinterpolation
  boundary/cf_interface/fortran/quadcfinterpolation
  boundary/physical_boundary/fortran/cartphysbdryop

  coarsen_ops/fortran/cubiccoarsen
  coarsen_ops/fortran/rt0coarsen
  lagrangian/fortran/lagrangian_interaction

  math/fortran/curl
  math/fortran/div
  math/fortran/flux
  math/fortran/grad
  math/fortran/graddetect
  math/fortran/interp
  math/fortran/laplace
  math/fortran/miscmath
  math/fortran/rot
  math/fortran/strain
  math/fortran/vclaplace

  refine_ops/fortran/divpreservingrefine
  refine_ops/fortran/cart_side_refine

  solvers/impls/fortran/patchsmoothers
  )

FOREACH(_d ${IBAMR_DIMENSIONS})
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    IBAMR_PROCESS_M4(NAME ${_fortran_root_name}${_d}d.f)
  ENDFOREACH()

  SET("FORTRAN_GENERATED_SRC${_d}D" "")
  FOREACH(_fortran_root_name ${FORTRAN_SRC_BASE})
    LIST(APPEND "FORTRAN_GENERATED_SRC${_d}D" "${CMAKE_CURRENT_BINARY_DIR}/${_fortran_root_name}${_d}d.f")
  ENDFOREACH()
ENDFOREACH()

#
# set up C++:
#
SET(CXX_SRC
  # boundary
  boundary/HierarchyGhostCellInterpolation.cpp
  boundary/cf_interface/CartCellDoubleLinearCFInterpolation.cpp
  boundary/cf_interface/CartSideDoubleQuadraticCFInterpolation.cpp
  boundary/cf_interface/CartCellDoubleQuadraticCFInterpolation.cpp
  boundary/physical_boundary/CartCellRobinPhysBdryOp.cpp
  boundary/physical_boundary/StaggeredPhysicalBoundaryHelper.cpp
  boundary/physical_boundary/CartExtrapPhysBdryOp.cpp
  boundary/physical_boundary/CartSideRobinPhysBdryOp.cpp
  boundary/physical_boundary/PhysicalBoundaryUtilities.cpp
  boundary/physical_boundary/muParserRobinBcCoefs.cpp
  boundary/physical_boundary/ExtendedRobinBcCoefStrategy.cpp
  boundary/physical_boundary/RobinPhysBdryPatchStrategy.cpp

  # coarsen ops
  coarsen_ops/CartCellDoubleCubicCoarsen.cpp
  coarsen_ops/CartSideDoubleRT0Coarsen.cpp
  coarsen_ops/LMarkerCoarsen.cpp
  coarsen_ops/CartSideDoubleCubicCoarsen.cpp

  # lagrangian
  lagrangian/LDataManager.cpp
  lagrangian/LNode.cpp
  lagrangian/LTransaction.cpp
  lagrangian/LEInteractor.cpp
  lagrangian/LNodeIndex.cpp
  lagrangian/LIndexSetData.cpp
  lagrangian/LSet.cpp
  lagrangian/LIndexSetDataFactory.cpp
  lagrangian/LSetData.cpp
  lagrangian/LIndexSetVariable.cpp
  lagrangian/LSetDataFactory.cpp
  lagrangian/LInitStrategy.cpp
  lagrangian/LSetDataIterator.cpp
  lagrangian/LMarker.cpp
  lagrangian/LSetVariable.cpp
  lagrangian/LData.cpp
  lagrangian/LMesh.cpp
  lagrangian/LSiloDataWriter.cpp

  # math
  math/PETScVecUtilities.cpp
  math/HierarchyMathOps.cpp
  math/PoissonUtilities.cpp
  math/SAMRAIGhostDataAccumulator.cpp
  math/PETScMatUtilities.cpp
  math/PatchMathOps.cpp

  # refine ops
  refine_ops/CartCellDoubleQuadraticRefine.cpp
  refine_ops/LMarkerRefine.cpp
  refine_ops/CartSideDoubleSpecializedLinearRefine.cpp
  refine_ops/CartCellDoubleBoundsPreservingConservativeLinearRefine.cpp
  refine_ops/CartSideDoubleDivPreservingRefine.cpp
  refine_ops/CartSideDoubleRT0Refine.cpp

  # solvers
  solvers/solver_utilities.cpp
  solvers/interfaces/FACPreconditionerStrategy.cpp
  solvers/interfaces/LinearOperator.cpp
  solvers/interfaces/JacobianOperator.cpp
  solvers/interfaces/GeneralOperator.cpp
  solvers/interfaces/GeneralSolver.cpp
  solvers/interfaces/KrylovLinearSolver.cpp
  solvers/interfaces/NewtonKrylovSolver.cpp
  solvers/interfaces/LinearSolver.cpp
  solvers/wrappers/PETScSNESFunctionGOWrapper.cpp
  solvers/wrappers/PETScSAMRAIVectorReal.cpp
  solvers/wrappers/PETScMatLOWrapper.cpp
  solvers/wrappers/PETScSNESJacobianJOWrapper.cpp
  solvers/wrappers/PETScPCLSWrapper.cpp
  solvers/impls/CCPoissonPETScLevelSolver.cpp
  solvers/impls/PoissonSolver.cpp
  solvers/impls/NewtonKrylovSolverManager.cpp
  solvers/impls/VCSCViscousPETScLevelSolver.cpp
  solvers/impls/PETScLevelSolver.cpp
  solvers/impls/CCPoissonSolverManager.cpp
  solvers/impls/PoissonFACPreconditionerStrategy.cpp
  solvers/impls/BGaussSeidelPreconditioner.cpp
  solvers/impls/BJacobiPreconditioner.cpp
  solvers/impls/SCPoissonPointRelaxationFACOperator.cpp
  solvers/impls/KrylovLinearSolverPoissonSolverInterface.cpp
  solvers/impls/VCSCViscousOperator.cpp
  solvers/impls/FACPreconditioner.cpp
  solvers/impls/PETScMFFDJacobianOperator.cpp
  solvers/impls/VCSCViscousOpPointRelaxationFACOperator.cpp
  solvers/impls/CCLaplaceOperator.cpp
  solvers/impls/PETScNewtonKrylovSolver.cpp
  solvers/impls/CCPoissonBoxRelaxationFACOperator.cpp
  solvers/impls/CCPoissonHypreLevelSolver.cpp
  solvers/impls/PoissonFACPreconditioner.cpp
  solvers/impls/CCPoissonLevelRelaxationFACOperator.cpp
  solvers/impls/PETScKrylovLinearSolver.cpp
  solvers/impls/KrylovLinearSolverManager.cpp
  solvers/impls/SCPoissonPETScLevelSolver.cpp
  solvers/impls/SCPoissonSolverManager.cpp
  solvers/impls/LaplaceOperator.cpp
  solvers/impls/SCLaplaceOperator.cpp
  solvers/impls/SCPoissonHypreLevelSolver.cpp
  solvers/impls/PETScKrylovPoissonSolver.cpp
  solvers/impls/CCPoissonPointRelaxationFACOperator.cpp

  # utilities
  utilities/RefinePatchStrategySet.cpp
  utilities/NodeSynchCopyFillPattern.cpp
  utilities/ParallelEdgeMap.cpp
  utilities/DebuggingUtilities.cpp
  utilities/CellNoCornersFillPattern.cpp
  utilities/CoarsenPatchStrategySet.cpp
  utilities/SideDataSynchronization.cpp
  utilities/StandardTagAndInitStrategySet.cpp
  utilities/IndexUtilities.cpp
  utilities/ParallelSet.cpp
  utilities/FaceDataSynchronization.cpp
  utilities/HierarchyIntegrator.cpp
  utilities/MergingLoadBalancer.cpp
  utilities/CopyToRootSchedule.cpp
  utilities/AppInitializer.cpp
  utilities/IBTKInit.cpp
  utilities/SAMRAIDataCache.cpp
  utilities/FixedSizedStream.cpp
  utilities/muParserCartGridFunction.cpp
  utilities/IBTK_MPI.cpp
  utilities/ibtk_utilities.cpp
  utilities/FaceSynchCopyFillPattern.cpp
  utilities/CartGridFunctionSet.cpp
  utilities/EdgeDataSynchronization.cpp
  utilities/NodeDataSynchronization.cpp
  utilities/ParallelMap.cpp
  utilities/SideNoCornersFillPattern.cpp
  utilities/box_utilities.cpp
  utilities/Streamable.cpp
  utilities/CopyToRootTransaction.cpp
  utilities/SideSynchCopyFillPattern.cpp
  utilities/CartGridFunction.cpp
  utilities/NormOps.cpp
  utilities/EdgeSynchCopyFillPattern.cpp
  utilities/SecondaryHierarchy.cpp
  utilities/StreamableManager.cpp
  utilities/LMarkerUtilities.cpp
  utilities/PartitioningBox.cpp
  )

IF(IBAMR_HAVE_LIBMESH)
  LIST(APPEND CXX_SRC
    # lagrangian
    lagrangian/BoxPartitioner.cpp
    lagrangian/FEDataInterpolation.cpp
    lagrangian/FEDataManager.cpp
    lagrangian/FEMapping.cpp
    lagrangian/FEProjector.cpp
    lagrangian/FEValues.cpp
    lagrangian/FischerGuess.cpp
    lagrangian/StableCentroidPartitioner.cpp

    # utilities
    utilities/LibMeshSystemVectors.cpp
    utilities/LibMeshSystemIBVectors.cpp
    utilities/libmesh_utilities.cpp
    )
ENDIF()

# TARGET_SOURCES needs to use private - if we used public then dependencies on
# this target will also compile the source files
TARGET_SOURCES(IBTK2d PRIVATE ${FORTRAN_GENERATED_SRC2D} ${CXX_SRC})
TARGET_SOURCES(IBTK3d PRIVATE ${FORTRAN_GENERATED_SRC3D} ${CXX_SRC})

# Include path to ibtk/config.h - note that IBTKHeaders covers the normal header
# locations both for building and installation
TARGET_INCLUDE_DIRECTORIES(IBTK2d PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/ibtk/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>)
TARGET_INCLUDE_DIRECTORIES(IBTK3d PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/ibtk/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/>)

TARGET_LINK_LIBRARIES(IBTK2d PUBLIC IBTKHeaders)
TARGET_LINK_LIBRARIES(IBTK3d PUBLIC IBTKHeaders)

IBAMR_SETUP_TARGET_LIBRARY(IBTK2d)
IBAMR_SETUP_TARGET_LIBRARY(IBTK3d)

INSTALL(TARGETS IBTK2d EXPORT IBAMRTargets COMPONENT library)
INSTALL(TARGETS IBTK3d EXPORT IBAMRTargets COMPONENT library)
