// Filename: INSVCStaggeredConservativeConvectiveOperator.h
// Created on 01 April 2018 by Nishant Nangia and Amneet Bhalla
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
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

#ifndef included_IBAMR_INSVCStaggeredConservativeConvectiveOperator
#define included_IBAMR_INSVCStaggeredConservativeConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "SideVariable.h"
#include "boost/array.hpp"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSVCStaggeredConservativeConvectiveOperator is a concrete
 * ConvectiveOperator that implements an upwind convective differencing operator.
 *
 * Class INSVCStaggeredConservativeConvectiveOperator computes the convective derivative of
 * a side-centered velocity field using various limiters described by Patel and Natarajan.
 *
 * A side-centered density update is provided by this class, which is used in the conservative discretization
 * of the Stokes operator. There is an optional mass density source term that can be set to easily check order
 * of accuracy via manufactured solutions.
 *
 * References
 * Patel, JK. and Natarajan, G., <A HREF="https://www.sciencedirect.com/science/article/pii/S0045793014004009">
 * A generic framework for design of interface capturing schemes for multi-fluid flows</A>
 *
 * \note This class is specialized in that it computes a conservative discretization of the form
 * \f$N = \nabla \cdot (u \rho u)\f$, where the density \f$\rho\f$ can vary in space and time.
 * This operator is to be used in conjuction with the conservative form of the variable coefficient
 * Navier-Stokes equations, which will produce better results for high density ratio flows.
 *
 * \see INSVCStaggeredHierarchyIntegrator
 */
class INSVCStaggeredConservativeConvectiveOperator : public ConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSVCStaggeredConservativeConvectiveOperator(const std::string& object_name,
                                                 SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                 ConvectiveDifferencingType difference_form,
                                                 std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> bc_coefs);

    /*!
     * \brief Destructor.
     */
    ~INSVCStaggeredConservativeConvectiveOperator();

    /*!
     * \brief Static function to construct an INSVCStaggeredConservativeConvectiveOperator.
     */
    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      ConvectiveDifferencingType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    {
        return new INSVCStaggeredConservativeConvectiveOperator(object_name, input_db, difference_form, bc_coefs);
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
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);

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

    /*!
     * \brief Set the current side-centered density patch data index.
     */
    void setSideCenteredDensityPatchDataIndex(int rho_sc_idx);

    /*
     * \brief Set the boundary condition object for the side-centered density.
     */
    void setSideCenteredDensityBoundaryConditions(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& rho_sc_bc_coefs);

    /*!
     * \brief Get the newly constructed side-centered density patch data index.
     *
     * \note This data is produced as a part of the apply() routine and should be used
     * in the linear operator for the INSVC solver.
     */
    int getUpdatedSideCenteredDensityPatchDataIndex();

    /*!
     * \brief Set an optional mass density source term.
     */
    void setMassDensitySourceTerm(const SAMRAI::tbox::Pointer<IBTK::CartGridFunction> S_fcn);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredConservativeConvectiveOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredConservativeConvectiveOperator(const INSVCStaggeredConservativeConvectiveOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredConservativeConvectiveOperator&
    operator=(const INSVCStaggeredConservativeConvectiveOperator& that) = delete;

    /*!
     * \brief Compute the advection velocity using simple averages
     */
    void computeAdvectionVelocity(
        boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > U_data,
        const SAMRAI::hier::IntVector<NDIM>& patch_lower,
        const SAMRAI::hier::IntVector<NDIM>& patch_upper,
        const boost::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes);

    /*!
     * \brief Compute the interpolation of a quantity Q onto Q_half, faces of the velocity DOF centered control volumes
     */
    void interpolateSideQuantity(
        boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> Q_half_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > Q_data,
        const SAMRAI::hier::IntVector<NDIM>& patch_lower,
        const SAMRAI::hier::IntVector<NDIM>& patch_upper,
        const boost::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const LimiterType& convective_limiter);

    /*!
     * \brief Compute div[rho_half*u_half*u_adv]
     */
    void computeConvectiveDerivative(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > N_data,
        boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> P_half_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> R_half_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_half_data,
        const boost::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const double* const dx);

    /*!
     * \brief Compute the density update rho = a0*rho^0 + a1*rho^1 + a2*dt*(-div[u_adv*rho_half]) + a2*dt*S
     */
    void computeDensityUpdate(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > R_data,
        const double& a0,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > R0_data,
        const double& a1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > R1_data,
        const double& a2,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> R_half_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > S_data,
        const boost::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const double& dt,
        const double* const dx);

    // Boundary condition helper object.
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;

    // Cached communications operators.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    std::string d_bdry_extrap_type;
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    // Whether or not the current density field has been set.
    bool d_rho_is_set;

    // Number of RK steps to take.
    int d_num_steps;

    // Boundary condition object for side-centered density field.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_rho_sc_bc_coefs;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_var;
    int d_U_scratch_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_rho_sc_var;
    int d_rho_sc_current_idx, d_rho_sc_scratch_idx, d_rho_sc_new_idx;

    // Hierarchy operation obect.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;

    // The limiter type for interpolation onto faces.
    LimiterType d_velocity_convective_limiter;
    LimiterType d_density_convective_limiter;

    // The required number of ghost cells for the chosen interpolation
    int d_velocity_limiter_gcw, d_density_limiter_gcw;

    // Variable to indicate the density update time-stepping type.
    TimeSteppingType d_density_time_stepping_type;

    // Source term variable and function for the mass density update
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_S_var;
    int d_S_scratch_idx;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_S_fcn;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSVCStaggeredConservativeConvectiveOperator
