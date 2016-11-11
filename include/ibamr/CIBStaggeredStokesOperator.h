// Filename: CIBStaggeredStokesOperator.h
// Created on 31 Oct 2014 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

#ifndef included_IBAMR_CIBStaggeredStokesOperator
#define included_IBAMR_CIBStaggeredStokesOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "ibamr/StaggeredStokesOperator.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

namespace IBAMR
{
class CIBStrategy;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CIBStaggeredStokesOperator is a concrete IBTK::LinearOperator which
 * implements a staggered-grid (MAC) discretization of the incompressible Stokes
 * operator while maintaining the constraint of rigidity for the immersed structures.
 *
 * This class is intended to be used with an iterative (Krylov or Newton-Krylov)
 * incompressible flow solver.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class CIBStaggeredStokesOperator : public IBAMR::StaggeredStokesOperator
{
    //////////////////////////////////////////////////////////////////////////////
public:
    /*!
     * \brief Class constructor.
     */
    CIBStaggeredStokesOperator(const std::string& object_name,
                               SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> cib_strategy,
                               bool homogeneous_bc = true);

    /*!
     * \brief Destructor.
     */
    ~CIBStaggeredStokesOperator();

    //\{ // Operator functionality of IBAMR::StaggeredStokesOperator class.
    /*!
     * \brief Compute hierarchy dependent data required for computing y=Ax.
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);
    //\}

    //\{ // Additional functionality of CIBStaggeredStokesOperator.

    /*!
     * \name Linear operator functionality.
     */
    using LinearOperator::apply;
    virtual void apply(Vec x, Vec y);

    // Set scaling factors for various operators to improve the condition number
    // of the system.

    /*!
     * \brief Set scale factor for interp operator.
     */
    void setInterpScaleFactor(const double beta);

    /*!
     * \brief Set scale factor for spread operator.
     */
    void setSpreadScaleFactor(const double gamma);

    /*!
     * \brief Set scale factor for regularizing mobility matrix.
     */
    void setRegularizeMobilityFactor(const double delta);

    /*!
     * \brief Set if the mean of the Lagrangian force is to be subtracted
     * from the Eulerian force variable.
     *
     * \note This operation is needed for certain situations like Stokes flow
     * with periodic BCs.
     */
    void setNormalizeSpreadForce(const bool normalize_force);

    /*
     * Set y := y - A*0, i.e., shift the right-hand-side vector to account for
     * inhomogeneous boundary conditions.
     */
    using LinearOperator::modifyRhsForBcs;
    virtual void modifyRhsForBcs(Vec y);

    /*!
     * \brief Impose boudary conditions in the solution vector.
     */
    using LinearOperator::imposeSolBcs;
    virtual void imposeSolBcs(Vec x);
    //\}

    //////////////////////////////////////////////////////////////////////////////
protected:
    // Pointer to a constraint based rigid IB Method.
    SAMRAI::tbox::Pointer<IBAMR::CIBStrategy> d_cib_strategy;

    // Scaling factors for various operators.
    double d_scale_interp, d_scale_spread, d_reg_mob_factor;
    bool d_normalize_spread_force;
    //////////////////////////////////////////////////////////////////////////////

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CIBStaggeredStokesOperator(const CIBStaggeredStokesOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CIBStaggeredStokesOperator& operator=(const CIBStaggeredStokesOperator& that);
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_CIBStaggeredStokesOperator
