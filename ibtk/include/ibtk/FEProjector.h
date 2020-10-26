// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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
 */
class FEProjector
{
public:
    /// Constructor.
    FEProjector(libMesh::EquationSystems* equation_systems, bool enable_logging = true);

    /// Alternative constructor that takes in a shared pointer to an FEData object.
    FEProjector(std::shared_ptr<FEData> fe_data, bool enable_logging = true);

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
     * For more information on the algorithm used to generate the (nearly)
     * diagonal matrix see 'A note on mass lumping and related processes in the
     * finite element method', Hinton, 1976.
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
    buildStabilizedL2ProjectionSolver(const std::string& system_name, double epsilon, const std::map<unsigned int, double>* epsilon_map_ptr = nullptr);

    /*!
     * \return Pointer to vector representation of diagonal L2 mass matrix.
     * Unlike buildLumpedL2ProjectionSolver this matrix is always diagonal and
     * is always assembled without applying constraints.
     *
     * For more information on the algorithm used to generate the diagonal
     * matrix see 'A note on mass lumping and related processes in the finite
     * element method', Hinton, 1976.
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

private:
    /*!
     * Whether or not to log data to the screen: see
     * FEProjector::setLoggingEnabled() and
     * FEProjector::getLoggingEnabled().
     */
    bool d_enable_logging = false;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FEProjector
