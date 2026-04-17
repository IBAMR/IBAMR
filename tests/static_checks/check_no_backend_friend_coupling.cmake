set(src_root "${CMAKE_SOURCE_DIR}")
if(NOT EXISTS "${src_root}/ibtk/include/ibtk/PETScLevelSolver.h")
  get_filename_component(src_root "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
endif()

set(petsc_level_solver_header "${src_root}/ibtk/include/ibtk/PETScLevelSolver.h")
set(stokes_level_solver_header "${src_root}/include/ibamr/StaggeredStokesPETScLevelSolver.h")

set(files_to_scan
    "${src_root}/ibtk/include/ibtk/private/PETScLevelSolverPetscShellBackend.h"
    "${src_root}/ibtk/include/ibtk/private/PETScLevelSolverBlasLapackShellBackend.h"
    "${src_root}/ibtk/include/ibtk/private/PETScLevelSolverEigenShellBackendBase.h"
    "${src_root}/ibtk/src/solvers/impls/PETScLevelSolverPetscShellBackend.cpp"
    "${src_root}/ibtk/src/solvers/impls/PETScLevelSolverBlasLapackShellBackend.cpp"
    "${src_root}/ibtk/src/solvers/impls/PETScLevelSolverEigenPseudoinverseShellBackend.cpp"
    "${src_root}/ibtk/src/solvers/impls/PETScLevelSolverEigenReferenceShellBackend.cpp"
    "${src_root}/src/navier_stokes/StaggeredStokesEigenSchurComplementShellBackend.cpp")

file(READ "${petsc_level_solver_header}" petsc_header_text)
if(petsc_header_text MATCHES "friend class[ \t\r\n]+.*ShellBackend")
  message(FATAL_ERROR
    "PETScLevelSolver backend friend coupling detected in ${petsc_level_solver_header}")
endif()

file(READ "${stokes_level_solver_header}" stokes_header_text)
if(stokes_header_text MATCHES "friend class[ \t\r\n]+StaggeredStokesEigenSchurComplementShellBackend")
  message(FATAL_ERROR
    "StaggeredStokesPETScLevelSolver backend friend coupling detected in ${stokes_level_solver_header}")
endif()

foreach(scan_file IN LISTS files_to_scan)
  file(READ "${scan_file}" scan_text)
  if(scan_text MATCHES "d_solver\\.d_")
    message(FATAL_ERROR
      "Backend implementation uses PETScLevelSolver private members directly: ${scan_file}")
  endif()
endforeach()
