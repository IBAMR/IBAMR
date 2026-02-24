// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_StaggeredStokesPETScLevelSolver
#define included_IBAMR_StaggeredStokesPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StaggeredStokesSolver.h"

#include "ibtk/PETScLevelSolver.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAICellVariable.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRefineSchedule.h"
#include "SAMRAISAMRAIVectorReal.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIVariableContext.h"

#include "petscvec.h"

#include <set>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesPETScLevelSolver is a concrete PETScLevelSolver
 * for a staggered-grid (MAC) discretization of the incompressible Stokes
 * equations.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class StaggeredStokesPETScLevelSolver : public IBTK::PETScLevelSolver, public StaggeredStokesSolver
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesPETScLevelSolver(const std::string& object_name,
                                    SAMRAIPointer<SAMRAIDatabase> input_db,
                                    const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesPETScLevelSolver();

    /*!
     * \brief Static function to construct a StaggeredStokesPETScLevelSolver.
     */
    static SAMRAIPointer<StaggeredStokesSolver> allocate_solver(const std::string& object_name,
                                                                SAMRAIPointer<SAMRAIDatabase> input_db,
                                                                const std::string& default_options_prefix)
    {
        return new StaggeredStokesPETScLevelSolver(object_name, input_db, default_options_prefix);
    } // allocate_solver

protected:
    /*!
     * \brief Generate IS/subdomains for Schwartz type preconditioners.
     */
    void generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                               std::vector<std::set<int>>& nonoverlap_is) override;

    /*!
     * \brief Generate IS/subdomains for fieldsplit type preconditioners.
     */
    void generateFieldSplitSubdomains(std::vector<std::string>& field_names,
                                      std::vector<std::set<int>>& field_is) override;

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     */
    void initializeSolverStateSpecialized(const SAMRAISAMRAIVectorReal<double>& x,
                                          const SAMRAISAMRAIVectorReal<double>& b) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverStateSpecialized().
     */
    void deallocateSolverStateSpecialized() override;

    /*!
     * \brief Copy a generic vector to the PETSc representation.
     */
    void copyToPETScVec(Vec& petsc_x, SAMRAISAMRAIVectorReal<double>& x) override;

    /*!
     * \brief Copy a generic vector from the PETSc representation.
     */
    void copyFromPETScVec(Vec& petsc_x, SAMRAISAMRAIVectorReal<double>& x) override;

    /*!
     * \brief Copy solution and right-hand-side data to the PETSc
     * representation, including any modifications to account for boundary
     * conditions.
     */
    void setupKSPVecs(Vec& petsc_x,
                      Vec& petsc_b,
                      SAMRAISAMRAIVectorReal<double>& x,
                      SAMRAISAMRAIVectorReal<double>& b) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesPETScLevelSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesPETScLevelSolver(const StaggeredStokesPETScLevelSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesPETScLevelSolver& operator=(const StaggeredStokesPETScLevelSolver& that) = delete;

    /*!
     * \name PETSc objects.
     */
    //\{

    SAMRAIPointer<SAMRAIVariableContext> d_context;
    std::vector<int> d_num_dofs_per_proc;
    int d_u_dof_index_idx = IBTK::invalid_index, d_p_dof_index_idx = IBTK::invalid_index;
    int d_u_nullspace_idx = IBTK::invalid_index, d_p_nullspace_idx = IBTK::invalid_index;
    SAMRAIPointer<SAMRAISideVariable<int>> d_u_dof_index_var;
    SAMRAIPointer<SAMRAISideVariable<double>> d_u_nullspace_var;
    SAMRAIPointer<SAMRAICellVariable<int>> d_p_dof_index_var;
    SAMRAIPointer<SAMRAICellVariable<double>> d_p_nullspace_var;
    SAMRAIPointer<SAMRAIRefineSchedule> d_data_synch_sched, d_ghost_fill_sched;

    //\}
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesPETScLevelSolver
