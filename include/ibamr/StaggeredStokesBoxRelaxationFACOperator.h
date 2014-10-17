// Filename: StaggeredStokesBoxRelaxationFACOperator.h
// Created on 11 Jun 2010 by Boyce Griffith
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

#ifndef included_StaggeredStokesBoxRelaxationFACOperator
#define included_StaggeredStokesBoxRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>
#include <vector>

#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
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
 * \brief Class StaggeredStokesBoxRelaxationFACOperator is a concrete
 * StaggeredStokesFACPreconditionerStrategy implementing a box relaxation
 * (Vanka-type) smoother for use as a multigrid preconditioner.
*/
class StaggeredStokesBoxRelaxationFACOperator : public StaggeredStokesFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesBoxRelaxationFACOperator(const std::string& object_name,
                                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                            const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesBoxRelaxationFACOperator();

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

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
    StaggeredStokesBoxRelaxationFACOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesBoxRelaxationFACOperator(const StaggeredStokesBoxRelaxationFACOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesBoxRelaxationFACOperator& operator=(const StaggeredStokesBoxRelaxationFACOperator& that);

    /*
     * Box operator data.
     */
    std::vector<Mat> d_box_op;
    std::vector<Vec> d_box_e, d_box_r;
    std::vector<KSP> d_box_ksp;

    /*
     * Mappings from patch indices to patch operators.
     */
    std::vector<std::vector<boost::array<SAMRAI::hier::BoxList<NDIM>, NDIM> > > d_patch_side_bc_box_overlap;
    std::vector<std::vector<SAMRAI::hier::BoxList<NDIM> > > d_patch_cell_bc_box_overlap;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StaggeredStokesBoxRelaxationFACOperator
