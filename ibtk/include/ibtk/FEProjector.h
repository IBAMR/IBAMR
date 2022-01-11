// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_FEProjector
#define included_IBTK_FEProjector

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include <ibtk/FischerGuess.h>

#include <tbox/Pointer.h>
#include <tbox/Timer.h>

#include <libmesh/equation_systems.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/petsc_vector.h>

#include <map>
#include <string>

namespace IBTK
{
class FEData;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FEProjector coordinates data structures for projecting
 * fields in FE models.
 *
 * <h2>Parameters read from the input database</h2>
 * <ol>
 *   <li>num_fischer_vectors: Number of previous solution and RHS pairs stored for
 *   use in computing the initial guess to each linear system. Defaults to five.
 *   Using more vectors will require additional computational work when computing
 *   the initial guess (roughly N*N dot products, where N is the number of stored
 *   vector pairs) but will decrease the number of solver iterations. The default
 *   value lowers the number of solver iterations to, typically, no more than two
 *   or three, so its usually the right value.</li>
 * </ol>
 */
class FEProjector
{
public:
    /// Constructor.
    FEProjector(libMesh::EquationSystems* equation_systems,
                const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db);

    /// Alternative constructor that takes in a shared pointer to an FEData object.
    FEProjector(std::shared_ptr<FEData> fe_data, const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db);

    /// Deleted default constructor.
    FEProjector() = delete;

    /// Deleted copy constructor.
    FEProjector(const FEProjector& from) = delete;

    /// Deleted assignment operator.
    FEProjector& operator=(const FEProjector& that) = delete;

    /// Defaulted destructor.
    ~FEProjector() = default;

    /*!
     * \return Pointers to a linear solver and sparse matrix corresponding to a
     * L2 projection operator.
     */
    std::pair<libMesh::PetscLinearSolver<double>*, libMesh::PetscMatrix<double>*>
    buildL2ProjectionSolver(const std::string& system_name);

    /*!
     * \return Pointers to a linear solver and sparse matrix corresponding to a
     * L2 projection operator generated with mass lumping. This system will be
     * diagonal if the finite element field does not have any hanging nodes:
     * i.e., the difference between this function and buildDiagonalL2MassMatrix
     * is that this function assembles the mass matrix while imposing
     * constraints.
     *
     * This mass matrix is constructed using nodal quadrature.
     */
    std::pair<libMesh::PetscLinearSolver<double>*, libMesh::PetscMatrix<double>*>
    buildLumpedL2ProjectionSolver(const std::string& system_name);

    /*!
     * \return Pointers to a linear solver and sparse matrix corresponding to a
     * L2 projection operator with a local projection stabilization term.
     *
     * This local stabilization approach augments the standard projection with a
     * stabilization of the form (p, q) + epsilon (p - Pi p, q - Pi q) = (f, q).
     * Currently the only supported local projection is to piecewise constant
     * functions.
     *
     * Ref: https://epubs.siam.org/doi/abs/10.1137/S0036142905444482
     */
    std::pair<libMesh::PetscLinearSolver<double>*, libMesh::PetscMatrix<double>*>
    buildStabilizedL2ProjectionSolver(const std::string& system_name, double epsilon);

    /*!
     * \return Pointer to vector representation of diagonal L2 mass matrix.
     * Unlike buildLumpedL2ProjectionSolver this matrix is always diagonal and
     * is always assembled without applying constraints.
     *
     * This mass matrix is constructed using nodal quadrature.
     */
    libMesh::PetscVector<double>* buildDiagonalL2MassMatrix(const std::string& system_name);

    /*!
     * \brief Set U to be the L2 projection of F.
     */
    bool computeL2Projection(libMesh::PetscVector<double>& U,
                             libMesh::PetscVector<double>& F,
                             const std::string& system_name,
                             bool consistent_mass_matrix = true,
                             bool close_U = true,
                             bool close_F = true,
                             double tol = 1.0e-6,
                             unsigned int max_its = 100);

    /*!
     * \brief Set U to be the L2 projection of F with a local projection
     * stabilization term.
     *
     * This local stabilization approach augments the standard projection with a
     * stabilization of the form (p, q) + epsilon (p - Pi p, q - Pi q) = (f, q).
     * Currently the only supported local projection is to piecewise constant
     * functions.
     *
     * Ref: https://epubs.siam.org/doi/abs/10.1137/S0036142905444482
     */
    bool computeStabilizedL2Projection(libMesh::PetscVector<double>& U,
                                       libMesh::PetscVector<double>& F,
                                       const std::string& system_name,
                                       double epsilon,
                                       bool close_U = true,
                                       bool close_F = true,
                                       double tol = 1.0e-6,
                                       unsigned int max_its = 100);

    /*!
     * \brief Enable or disable logging.
     */
    void setLoggingEnabled(bool enable_logging = true);

    /*!
     * \brief Determine whether logging is enabled or disabled.
     */
    bool getLoggingEnabled() const;

protected:
    /*!
     * FEData object that contains the libMesh data structures.
     *
     * @note multiple FEDataManager objects may use the same FEData object,
     * usually combined with different hierarchies.
     */
    std::shared_ptr<FEData> d_fe_data;

    /// Data structures for consistent mass matrices and related solvers.
    std::map<std::string, std::unique_ptr<libMesh::PetscMatrix<double> > > d_L2_proj_matrix;
    std::map<std::string, std::unique_ptr<libMesh::PetscLinearSolver<double> > > d_L2_proj_solver;

    /// Data structures for lumped mass matrices. These are computed in the same
    /// way as the normal mass matrix, except the quadrature rule used to
    /// compute matrix entries has it's points at the finite element nodes: i.e.,
    /// in the absence of constraints, the matrix is diagonal.
    ///
    /// Here we refer to the unconstrained matrix (which is always diagonal) as
    /// proj_matrix_diag and the constrained (should there be constraints) matrix
    /// as lumped_L2_proj_matrix.
    std::map<std::string, std::unique_ptr<libMesh::PetscMatrix<double> > > d_lumped_L2_proj_matrix;
    std::map<std::string, std::unique_ptr<libMesh::PetscLinearSolver<double> > > d_lumped_L2_proj_solver;
    std::map<std::string, std::unique_ptr<libMesh::PetscVector<double> > > d_diag_L2_proj_matrix;

    /// Data structures for consistent mass matrices and related solvers with local projection stabilization.
    std::map<std::string, std::map<double, std::unique_ptr<libMesh::PetscMatrix<double> > > > d_stab_L2_proj_matrix;
    std::map<std::string, std::map<double, std::unique_ptr<libMesh::PetscLinearSolver<double> > > >
        d_stab_L2_proj_solver;

    std::map<std::string, FischerGuess> d_initial_guesses;

private:
    /*!
     * Pointers for system-specific solver timers.
     *
     * SAMRAI stores timers in a single unsorted array. To keep timer lookups
     * quick we cache the pointers here.
     */
    std::map<std::string, SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> > d_linear_solve_system_timers;

    /*!
     * Whether or not to log data to the screen: see
     * FEProjector::setLoggingEnabled() and
     * FEProjector::getLoggingEnabled().
     */
    bool d_enable_logging = false;

    /*!
     * Number of vectors to use in the FischerGuess objects.
     */
    int d_num_fischer_vectors = 5;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_FEProjector
