// Filename: INSCollocatedPPMConvectiveOperator.h
// Created on 24 Aug 2011 by Boyce Griffith
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

#ifndef included_INSCollocatedPPMConvectiveOperator
#define included_INSCollocatedPPMConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "SAMRAI/tbox/Database.h"

namespace SAMRAI
{
namespace solv
{
template <class TYPE>
class SAMRAIVectorReal;

class RobinBcCoefStrategy;
} // namespace solv
namespace xfer
{

class CoarsenSchedule;

class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSCollocatedPPMConvectiveOperator is a concrete
 * ConvectiveOperator that implements a upwind convective differencing operator
 * based on the piecewise parabolic method (PPM).
 *
 * Class INSCollocatedPPMConvectiveOperator computes the convective derivative of
 * a cell-centered velocity field using the xsPPM7 method of Rider, Greenough,
 * and Kamm.
 *
 * \see INSCollocatedHierarchyIntegrator
 */
class INSCollocatedPPMConvectiveOperator : public ConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSCollocatedPPMConvectiveOperator(const std::string& object_name,
                                       boost::shared_ptr<SAMRAI::tbox::Database> input_db,
                                       ConvectiveDifferencingType difference_form,
                                       const std::vector<SAMRAI::solv::RobinBcCoefStrategy*>& bc_coefs);

    /*!
     * \brief Destructor.
     */
    ~INSCollocatedPPMConvectiveOperator();

    /*!
     * \brief Static function to construct an
     * INSCollocatedPPMConvectiveOperator.
     */
    static boost::shared_ptr<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      boost::shared_ptr<SAMRAI::tbox::Database> input_db,
                      ConvectiveDifferencingType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy*>& bc_coefs)
    {
        return boost::make_shared<INSCollocatedPPMConvectiveOperator>(object_name, input_db, difference_form, bc_coefs);
    } // allocate_operator

    /*!
     * \brief Compute the action of the convective operator.
     */
    void applyConvectiveOperator(int U_idx, int N_idx);

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for computing y=F[x] and
     * z=F[x]+y.
     *
     * The vector arguments for apply(), applyAdjoint(), etc, need not match
     * those for initializeOperatorState().  However, there must be a certain
     * degree of similarity, including
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the input and output vectors
     *
     * \note It is generally necessary to reinitialize the operator state when
     * the hierarchy configuration changes.
     *
     * It is safe to call initializeOperatorState() when the state is already
     * initialized.  In this case, the operator state is first deallocated and
     * then reinitialized.
     *
     * Conditions on arguments:
     * - input and output vectors must have same hierarchy
     * - input and output vectors must have same structure, depth, etc.
     *
     * Call deallocateOperatorState() to remove any data allocated by this
     * method.
     *
     * \see deallocateOperatorState
     *
     * \param in input vector
     * \param out output vector
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<double>& out);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * \note It is safe to call deallocateOperatorState() when the operator
     * state is already deallocated.
     *
     * \see initializeOperatorState
     */
    void deallocateOperatorState();

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSCollocatedPPMConvectiveOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSCollocatedPPMConvectiveOperator(const INSCollocatedPPMConvectiveOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSCollocatedPPMConvectiveOperator& operator=(const INSCollocatedPPMConvectiveOperator& that);

    // Data communication algorithms, operators, and schedules.
    boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> d_coarsen_alg;
    std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> > d_coarsen_scheds;
    boost::shared_ptr<SAMRAI::xfer::RefineAlgorithm> d_ghostfill_alg;
    boost::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> d_ghostfill_strategy;
    std::vector<boost::shared_ptr<SAMRAI::xfer::RefineSchedule> > d_ghostfill_scheds;
    std::string d_bdry_extrap_type;

    // Hierarchy configuration.
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    // Scratch data.
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > d_U_var;
    int d_U_scratch_idx;
    boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > d_u_extrap_var, d_u_flux_var;
    int d_u_extrap_idx, d_u_flux_idx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSCollocatedPPMConvectiveOperator
