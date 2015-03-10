// Filename: StaggeredStokesPETScLevelSolver.h
// Created on 08 Sep 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_StaggeredStokesPETScLevelSolver
#define included_StaggeredStokesPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibtk/PETScLevelSolver.h"
#include "petscvec.h"
#include "SAMRAI/tbox/Database.h"

namespace SAMRAI
{
namespace hier
{

class PatchLevel;
} // namespace hier
namespace solv
{
template <class TYPE>
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
                                    boost::shared_ptr<SAMRAI::tbox::Database> input_db,
                                    const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesPETScLevelSolver();

    /*!
     * \brief Static function to construct a StaggeredStokesPETScLevelSolver.
     */
    static boost::shared_ptr<StaggeredStokesSolver> allocate_solver(const std::string& object_name,
                                                                    boost::shared_ptr<SAMRAI::tbox::Database> input_db,
                                                                    const std::string& default_options_prefix)
    {
        return boost::make_shared<StaggeredStokesPETScLevelSolver>(object_name, input_db, default_options_prefix);
    } // allocate_solver

protected:
    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     */
    void initializeSolverStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<double>& x,
                                          const SAMRAI::solv::SAMRAIVectorReal<double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverStateSpecialized().
     */
    void deallocateSolverStateSpecialized();

    /*!
     * \brief Copy a generic vector to the PETSc representation.
     */
    void copyToPETScVec(Vec& petsc_x,
                        SAMRAI::solv::SAMRAIVectorReal<double>& x,
                        boost::shared_ptr<SAMRAI::hier::PatchLevel> patch_level);

    /*!
     * \brief Copy a generic vector from the PETSc representation.
     */
    void copyFromPETScVec(Vec& petsc_x,
                          SAMRAI::solv::SAMRAIVectorReal<double>& x,
                          boost::shared_ptr<SAMRAI::hier::PatchLevel> patch_level);

    /*!
     * \brief Copy solution and right-hand-side data to the PETSc
     * representation, including any modifications to account for boundary
     * conditions.
     */
    void setupKSPVecs(Vec& petsc_x,
                      Vec& petsc_b,
                      SAMRAI::solv::SAMRAIVectorReal<double>& x,
                      SAMRAI::solv::SAMRAIVectorReal<double>& b,
                      boost::shared_ptr<SAMRAI::hier::PatchLevel> patch_level);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesPETScLevelSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesPETScLevelSolver(const StaggeredStokesPETScLevelSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesPETScLevelSolver& operator=(const StaggeredStokesPETScLevelSolver& that);

    /*!
     * \name PETSc objects.
     */
    //\{

    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;
    std::vector<int> d_num_dofs_per_proc;
    int d_u_dof_index_idx, d_p_dof_index_idx;
    boost::shared_ptr<SAMRAI::pdat::SideVariable<int> > d_u_dof_index_var;
    boost::shared_ptr<SAMRAI::pdat::CellVariable<int> > d_p_dof_index_var;
    boost::shared_ptr<SAMRAI::xfer::RefineSchedule> d_data_synch_sched, d_ghost_fill_sched;

    //\}
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StaggeredStokesPETScLevelSolver
