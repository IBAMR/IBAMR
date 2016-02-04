// Filename: StaggeredStokesLevelRelaxationFACOperator.h
// Created on 16 Mar 2015 by Amneet Bhalla
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

#ifndef included_StaggeredStokesLevelRelaxationFACOperator
#define included_StaggeredStokesLevelRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>
#include <vector>

#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

namespace boost
{
template <class T, std::size_t N>
class array;
} // namespace boost
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoxList;
} // namespace hier
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesLevelRelaxationFACOperator is a concrete
 * StaggeredStokesFACPreconditionerStrategy implementing a level relaxation
 * smoother for use as a multigrid preconditioner.
*/
class StaggeredStokesLevelRelaxationFACOperator : public StaggeredStokesFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesLevelRelaxationFACOperator(const std::string& object_name,
                                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                              const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesLevelRelaxationFACOperator();

    /*!
     * \brief Static function to construct a StaggeredStokesFACPreconditioner with a
     * StaggeredStokesLevelRelaxationFACOperator FAC strategy.
     */
    static SAMRAI::tbox::Pointer<StaggeredStokesSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        SAMRAI::tbox::Pointer<StaggeredStokesFACPreconditionerStrategy> fac_operator =
            new StaggeredStokesLevelRelaxationFACOperator(
                object_name + "::StaggeredStokesLevelRelaxationFACOperator", input_db, default_options_prefix);
        return new StaggeredStokesFACPreconditioner(object_name, fac_operator, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Specify the smoother type.
     */
    void setSmootherType(const std::string& level_solver_type);

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are
     *being
     *performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are
     *being
     *performed
     */
    void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                     const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps);

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    void initializeOperatorStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                                            int coarsest_reset_ln,
                                            int finest_reset_ln);

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    void deallocateOperatorStateSpecialized(int coarsest_reset_ln, int finest_reset_ln);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesLevelRelaxationFACOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesLevelRelaxationFACOperator(const StaggeredStokesLevelRelaxationFACOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesLevelRelaxationFACOperator& operator=(const StaggeredStokesLevelRelaxationFACOperator& that);

    /*
     * Level solvers and solver parameters.
     */
    std::string d_level_solver_type, d_level_solver_default_options_prefix;
    double d_level_solver_abs_residual_tol, d_level_solver_rel_residual_tol;
    int d_level_solver_max_iterations;
    std::vector<SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesSolver> > d_level_solvers;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_level_solver_db;

    /*
     * Mappings from patch indices to patch operators.
     */
    std::vector<std::vector<boost::array<SAMRAI::hier::BoxList<NDIM>, NDIM> > > d_patch_side_bc_box_overlap;
    std::vector<std::vector<SAMRAI::hier::BoxList<NDIM> > > d_patch_cell_bc_box_overlap;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StaggeredStokesLevelRelaxationFACOperator
