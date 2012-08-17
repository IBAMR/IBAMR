// Filename: INSStaggeredPETScLevelSolver.h
// Created on 08 Sep 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
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

#ifndef included_INSStaggeredPETScLevelSolver
#define included_INSStaggeredPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDE
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/StokesSolver.h>

// IBTK INCLUDES
#include <ibtk/PETScLevelSolver.h>

// SAMRAI INCLUDES
#include <RefineSchedule.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredPETScLevelSolver is a concrete PETScLevelSolver for
 * a staggered-grid (MAC) discretization of the incompressible Stokes equations.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class INSStaggeredPETScLevelSolver
    : public IBTK::PETScLevelSolver,
      public StokesSolver
{
public:
    /*!
     * \brief Constructor.
     */
    INSStaggeredPETScLevelSolver(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredPETScLevelSolver();

    /*!
     * \brief Solver name string.
     */
    static const std::string SOLVER_TYPE_NAME;

    /*!
     * \brief Static function to construct a CCPoissonPETScLevelSolver.
     */
    static SAMRAI::tbox::Pointer<StokesSolver>
    allocate_solver(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
        {
            return new INSStaggeredPETScLevelSolver(object_name, input_db);
        }// allocate_solver

protected:
    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     */
    void
    initializeSolverStateSpecialized(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverStateSpecialized().
     */
    void
    deallocateSolverStateSpecialized();

    /*!
     * \brief Copy a generic vector to the PETSc representation.
     */
    void
    copyToPETScVec(
        Vec& petsc_x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);

    /*!
     * \brief Copy a generic vector from the PETSc representation.
     */
    void
    copyFromPETScVec(
        Vec& petsc_x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);

    /*!
     * \brief Copy solution and right-hand-side data to the PETSc
     * representation, including any modifications to account for boundary
     * conditions.
     */
    void
    setupKSPVecs(
        Vec& petsc_x,
        Vec& petsc_b,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredPETScLevelSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredPETScLevelSolver(
        const INSStaggeredPETScLevelSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredPETScLevelSolver&
    operator=(
        const INSStaggeredPETScLevelSolver& that);

    /*!
     * \name PETSc objects.
     */
    //\{

    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    std::vector<int> d_num_dofs_per_proc;
    int d_u_dof_index_idx, d_p_dof_index_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,int> > d_u_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,int> > d_p_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > d_data_synch_sched, d_ghost_fill_sched;

    //\}
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredPETScLevelSolver.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredPETScLevelSolver
