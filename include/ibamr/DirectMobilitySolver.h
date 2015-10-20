// Filename: DirectMobilitySolver.h
// Created on 20 Feb 2015 by Amneet Bhalla and Bakytzhan Kallemov
//
// Copyright (c) 2002-2015, Amneet Bhalla and Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_DirectMobilitySolver
#define included_DirectMobilitySolver

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <vector>
#include <map>
#include <string>

#include "petscvec.h"
#include "ibamr/ibamr_enums.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI
namespace IBAMR
{
class StokesSpecifications;
class CIBStrategy;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{

/*!
 * \brief Class DirectMobilitySolver solves the mobility sub-problem by
 * employing direct solvers for computing the inverse of mobility matrix.
 */
enum MobilitySubmatrixManager
{
    MM_EFFICIENT=0,
    MM_MEMORY_SAVE,
    MM_NONE,
    MM_UNKNOWN_TYPE = -1
};
enum  MobSubMatrixType
{
    BLOCK_DIAG=0,
    DENSE,
    CODE_CUSTOM,
    MSM_UNKNOWN_TYPE = -1
};
class DirectMobilitySolver : public SAMRAI::tbox::DescribedClass
{

    /////////////////////////////// PUBLIC //////////////////////////////////////
public:
    /*!
     * \brief The only constructor of this class.
     */
    DirectMobilitySolver(const std::string& object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy);

    /*!
     * \brief Destructor for this class.
     */
    ~DirectMobilitySolver();

    /*!
         * \brief Register a prototypical structure with a particular mobility
         * matrix identified by its name.
         *
         * \param mat_name Matrix handle.
         *
         * \param prototype_struct_id Prototypical structure id as defined while
         * registering structures with IBAMR::IBStrategy class.
         *
         * \param mat_type Matrix type to be used for dense mobility matrix.
         *
         * \param inv_type Inversion method to be used to invert the dense matrix.
         *
         * \param filename If the mobility matrix is to be read from an input file.
         *
         * \param scale Scale for improving the conditioning number of dense mobility
         * matrix. The matrix is scaled as \f$ [MM] = \alpha*[MM] + \beta*[I]. \f$
         *
         * \param managing_proc MPI processor that manages this mobility matrix.
         */
    void registerMobilityMat(const std::string& mat_name,
                             const unsigned prototype_struct_id,
                             MobilityMatrixType mat_type,
                             MobilityMatrixInverseType mob_inv_type,
                             MobilityMatrixInverseType body_inv_type,
                             const int managing_proc = 0,
                             const std::string& filename = "",
                             std::pair<double, double> scale = std::pair<double, double>(1.0, 0.0));

    /*!
     * \brief Register mutiple prototype structures with a particular matrix
     * identified by its name. A combined dense mobility matrix for all the
     * associated prototypical structures will be formed.
     *
     * \param mat_name Matrix handle.
     *
     * \param prototype_struct_ids Vector of prototypical structure ids as
     * defined while registering structures with IBAMR::IBStrategy class.
     *
     * \param mat_type Matrix type to be used for dense mobility matrix.
     *
     * \param inv_type Inversion method to be used to invert the dense matrix.
     *
     * \param scale Scale for improving the conditioning number of dense mobility
     * matrix. The matrix is scaled as \f$ [MM] = \alpha*[MM] + \beta*[I]. \f$
     *
     * \param managing_proc MPI processor that manages this mobility matrix.
     */
    void registerMobilityMat(const std::string& mat_name,
                             const std::vector<unsigned>& prototype_struct_ids,
                             MobilityMatrixType mat_type,
                             MobilityMatrixInverseType mob_inv_type,
                             MobilityMatrixInverseType body_inv_type,
                             const int managing_proc = 0,
                             const std::string& filename = "",
                             std::pair<double, double> scale = std::pair<double, double>(1.0, 0.0));

    /*!
     * \brief Register all structures that will be managed by this particular
     * dense mobility matrix.
     *
     * \note This function should be called after all the mobility matrices
     * have been registered with a particular processor.
     */
    void registerStructIDsWithMobilityMat(const std::string& mat_name,
                                          const std::vector<std::vector<unsigned> >& struct_ids);

    /*!
     * \brief Get access to friction (inverse of mobility) matrix.
     *
     * \param mat_name Matrix handle.
     *
     * \param fm Linear (fortran-style column major) array holding the matrix.
     *
     * \param size Column or row size of the square mobility matrix.
     *
     * \param managing_proc Rank of the processor managing the friction
     *  matrix.
     *
     *  \param inv_method Type of factorization used to invert the mobility matrix.
     */
    void getFrictionMat(const std::string& mat_name,
                        double** fm,
                        int* size,
                        int* managing_proc,
                        MobilityMatrixInverseType* inv_method);

    /*!
     * \brief Initialize the solver.
     */
    void initializeSolverState(Vec x, Vec b);

    /*!
     * \brief Deallocate the solver.
     */
    void deallocateSolverState();

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void setSolutionTime(const double solution_time);

    /*!
     * \brief Set the time interval of integration.
     */
    void setTimeInterval(const double current_time, const double new_time);

    /*!
     * \brief Set stokes specifications.
     */
    void setStokesSpecifications(const IBAMR::StokesSpecifications& stokes_spec);

    /*!
     * \brief Solves the mobility problem.
     *
     * \param x Vec storing the Lagrange multiplier
     *
     * \param b Vec storing the desired velocity
     *
     * \param skip_nonfree_parts Boolean indicating if the solution
     * is to be set only for free moving components.
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveSystem(Vec x, Vec b, const bool skip_nonfree_parts=false);

    /*!
     * \brief Solves the body mobility problem.
     *
     * \param x Vec storing the body velocity
     *
     * \param b Vec storing the applied force
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveBodySystem(Vec x, Vec b);

    /*!
     * \brief Return the ids of the structures associated with the dense
         * mobility matrix formation.
         *
         * \param mat_name Matrix handle.
         *
         * \return Vector of structure ids used to form the mobility matrix.
     */
    const std::vector<unsigned>& getPrototypeStructIDs(const std::string& mat_name);

    /*!
     * \brief Return the ids of the structure associated with the dense
     * mobility matrix.
     *
     * \param mat_name Matrix handle.
     *
     * \return Vector of structure ids associated with the mobility matrix.
     */
    const std::vector<std::vector<unsigned> >& getStructIDs(const std::string& mat_name);

    /////////////////////////////// PRIVATE //////////////////////////////////////
private:
    /*!
     * \brief Get input options.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief recreate mobility matrix.
     */
    void createMobilityMatrix();

    /*!
         * \brief Compute the inverse of mobility matrix using direct solvers.
         */
    void generateFrictionMatrix();

    /*!
         * \brief Compute the a body mobility matrix using direct solvers.
         */
    void generateBodyMobilityMatrix();

    /*!
     * \brief Compute solution and store in the rhs vector.
     */
    void computeSolution(const std::string& mat_name, double* managed_mat, const int mat_size, const MobilityMatrixInverseType inv_type, double* rhs, const bool BodyMobility=false, const int nrhs=1);

    /*!
     * \brief provide a matrix decomposition.
     */
    void decomposeMatrix(const std::string& mat_name, double* managed_mat, const int mat_size, const MobilityMatrixInverseType inv_type, const bool BodyMobility=false);

    /*!
     * \brief computes a matrix inverse.
     */
    void getInverseMatrix(const std::string& mat_name, double* managed_mat, double* inverse_mat, const int mat_size, const MobilityMatrixInverseType inv_type, const bool BodyMobility=false);

    /*!
     * \brief defines a manager for global mobility matrix.
     */
    void runMobilityMatManager();

    // Solver stuff
    std::string d_object_name;
    bool d_is_initialized;
    double d_solution_time, d_current_time, d_new_time;
    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_patch_hierarchy;

    // Structure(s) stuff.

    std::map<std::string, double*> d_managed_mat_map;
    std::map<std::string, double*> d_body_mob_mat_map;
    std::map<std::string, int> d_managed_mat_proc_map;
    std::map<std::string, std::vector<unsigned> > d_managed_mat_prototype_id_map;
    std::map<std::string, std::vector<std::vector<unsigned> > > d_managed_mat_actual_id_map;
    std::map<std::string, MobilityMatrixType> d_managed_mat_type_map;
    std::map<std::string, MobilityMatrixInverseType> d_managed_mat_inv_type_map;
    std::map<std::string, MobilityMatrixInverseType> d_body_mat_inv_type_map;
    std::map<std::string, unsigned int> d_managed_mat_nodes_map;
    std::map<std::string, std::pair<double, double> > d_managed_mat_scale_map;
    std::map<std::string, std::string> d_managed_mat_filename_map;
    std::map<unsigned, std::string > d_actual_id_managed_mat_map;

    // System physical parameters.
    double d_mu;
    double d_rho;

    // Parameters used in this class.
    std::string d_kernel_name;
    double d_f_periodic_corr;
    bool d_recompute_mob_mat;
    double d_svd_inv_tol, d_svd_inv_eps, d_body_svd_inv_tol, d_body_svd_inv_eps;
    std::map<std::string, int*> d_ipiv; // for LAPACK LU calls
    std::map<std::string, int*> d_body_ipiv; // for LAPACK LU calls
    bool d_mob_mat_register_completed;
    MobSubMatrixType d_submat_type;
    MobilitySubmatrixManager d_manager_type;
    MobilityMatrixType d_mob_mat_type;
    MobilityMatrixInverseType d_mob_mat_inv_type, d_body_mat_inv_type;
    std::pair<double, double> d_mob_scale;
    std::vector<std::string> d_stored_mob_mat_filenames;

}; // DirectMobilitySolver

} // IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_DirectMobilitySolver
