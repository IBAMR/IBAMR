// Filename: IntermediateVelocityRobinBcCoefs.C
// Last modified: <17.May.2007 13:17:20 griffith@box221.cims.nyu.edu>
// Created on 30 Sep 2006 by Boyce Griffith (boyce@trasnaform2.local)

#include "IntermediateVelocityRobinBcCoefs.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// STOOLS INCLUDES
#include <stools/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <ArrayDataBasicOps.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IntermediateVelocityRobinBcCoefs::IntermediateVelocityRobinBcCoefs(
    const int velocity_depth,
    const bool using_pressure_increment_form,
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const bc_coef)
    : d_velocity_depth(velocity_depth),
      d_using_pressure_increment_form(using_pressure_increment_form),
      d_bc_coef(bc_coef),
      d_correct_bc_coefs(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_rho(std::numeric_limits<double>::quiet_NaN()),
      d_P_var(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >(NULL)),
      d_P_src_idx(-1),
      d_P_dst_idx(-1),
      d_Phi_var(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >(NULL)),
      d_Phi_src_idx(-1),
      d_Phi_dst_idx(-1)
{
    // intentionally blank
    return;
}// IntermediateVelocityRobinBcCoefs

IntermediateVelocityRobinBcCoefs::~IntermediateVelocityRobinBcCoefs()
{
    // intentionally blank
    return;
}// ~IntermediateVelocityRobinBcCoefs

void
IntermediateVelocityRobinBcCoefs::correctBcCoefs(
    const bool correct_bc_coefs)
{
    d_correct_bc_coefs = correct_bc_coefs;
    return;
}// correctBcCoefs

void
IntermediateVelocityRobinBcCoefs::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time     =     new_time;
    return;
}// setCurrentTime

void
IntermediateVelocityRobinBcCoefs::setRho(
    const double rho)
{
    d_rho = rho;
    return;
}// setRho

void
IntermediateVelocityRobinBcCoefs::setPressureVariable(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_var,
    const int P_src_idx,
    const int P_dst_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!P_var.isNull());
    assert(P_src_idx != -1);
    assert(P_dst_idx != -1);
#endif
    d_P_var = P_var;
    d_P_src_idx = P_src_idx;
    d_P_dst_idx = P_dst_idx;
    return;
}// setPressureVariable

void
IntermediateVelocityRobinBcCoefs::setPhiVariable(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Phi_var,
    const int Phi_src_idx,
    const int Phi_dst_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!Phi_var.isNull());
    assert(Phi_src_idx != -1);
    assert(Phi_dst_idx != -1);
#endif
    d_Phi_var = Phi_var;
    d_Phi_src_idx = Phi_src_idx;
    d_Phi_dst_idx = Phi_dst_idx;
    return;
}// setPhiIndex

void
IntermediateVelocityRobinBcCoefs::getBcFillVars(
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > >& bc_fill_vars,
    std::vector<int>& bc_fill_src_idxs,
    std::vector<int>& bc_fill_dst_idxs) const
{
    if (d_velocity_depth == 0)
    {
        bc_fill_vars.resize(2);
        bc_fill_vars[0] =   d_P_var;
        bc_fill_vars[1] = d_Phi_var;

        bc_fill_src_idxs.resize(2);
        bc_fill_src_idxs[0] =   d_P_src_idx;
        bc_fill_src_idxs[1] = d_Phi_src_idx;

        bc_fill_dst_idxs.resize(2);
        bc_fill_dst_idxs[0] =   d_P_dst_idx;
        bc_fill_dst_idxs[1] = d_Phi_dst_idx;
    }
    else
    {
        bc_fill_vars.resize(0);
        bc_fill_src_idxs.resize(0);
        bc_fill_dst_idxs.resize(0);
    }
    return;
}// getBcFillVars

void
IntermediateVelocityRobinBcCoefs::setBcCoefs(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
#if USING_OLD_ROBIN_BC_INTERFACE
    TBOX_ERROR("IntermediateVelocityRobinBcCoefs::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#else
    d_bc_coef->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif
    if (d_correct_bc_coefs && fill_time > d_current_time)
    {
        computeCorrectedBcCoefs(acoef_data, bcoef_data, gcoef_data, patch, bdry_box);
    }
    return;
}// setBcCoefs

void
IntermediateVelocityRobinBcCoefs::setBcCoefs(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
#if USING_OLD_ROBIN_BC_INTERFACE
    d_bc_coef->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    TBOX_ERROR("IntermediateVelocityRobinBcCoefs::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#endif
    if (d_correct_bc_coefs && fill_time > d_current_time)
    {
        const SAMRAI::hier::Box<NDIM>& bc_coef_box = acoef_data->getBox();
        SAMRAI::math::ArrayDataBasicOps<NDIM,double> array_ops;
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > bcoef_data =
            new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
        array_ops.scale(*bcoef_data, -1.0, *acoef_data, bc_coef_box);
        array_ops.addScalar(*bcoef_data, *bcoef_data, 1.0, bc_coef_box);

        computeCorrectedBcCoefs(acoef_data, bcoef_data, gcoef_data, patch, bdry_box);
    }
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
IntermediateVelocityRobinBcCoefs::numberOfExtensionsFillable() const
{
    return d_bc_coef->numberOfExtensionsFillable();
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IntermediateVelocityRobinBcCoefs::computeCorrectedBcCoefs(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box) const
{
    // We only modify the inhomogeneous boundary data.
    if (gcoef_data.isNull()) return;

    // Get pointers to the patch data required to correct the boundary data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> >   P_data = patch.getPatchData(  d_P_dst_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Phi_data = patch.getPatchData(d_Phi_dst_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!  P_data.isNull());
    assert(!Phi_data.isNull());
#endif

    // Correct the boundary data.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double dt = d_new_time - d_current_time;

    const SAMRAI::hier::Box<NDIM> bc_coef_box = STOOLS::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_box);
    const int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis =  location_index / 2;
    const bool bdry_upper_side = (location_index % 2) != 0;

    for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const SAMRAI::hier::Index<NDIM>& i_s_bdry = b();
        const double& a = (*acoef_data)(i_s_bdry,0);
        const double& b = (*bcoef_data)(i_s_bdry,0);
        double& g = (*gcoef_data)(i_s_bdry,0);

        assert(b == 0.0);

        // i_s_bdry: side index located on physical boundary
        //
        // i_c_intr: cell index located adjacent to physical boundary in the
        // patch interior
        //
        // i_c_bdry: cell index located adjacent to physical boundary in the
        // patch exterior
        SAMRAI::pdat::CellIndex<NDIM> i_c_intr = i_s_bdry;
        SAMRAI::pdat::CellIndex<NDIM> i_c_bdry = i_c_intr;
        if (bdry_upper_side)
        {
            i_c_intr(bdry_normal_axis) -= 1;
            i_c_bdry(bdry_normal_axis) = i_c_intr(bdry_normal_axis)+1;
        }
        else
        {
            i_c_bdry(bdry_normal_axis) = i_c_intr(bdry_normal_axis)-1;
        }

        if (d_velocity_depth == bdry_normal_axis)
        {
            // Correct boundary conditions for the normal velocity.
            const double   dP_dn = ((*  P_data)(i_c_bdry)-(*  P_data)(i_c_intr))/dx[bdry_normal_axis];
            const double dPhi_dn = ((*Phi_data)(i_c_bdry)-(*Phi_data)(i_c_intr))/dx[bdry_normal_axis];

            if (d_using_pressure_increment_form)
            {
                if (bdry_upper_side)
                {
                    g -= a*(dt/d_rho)*(0.0  -dPhi_dn);
                }
                else
                {
                    g += a*(dt/d_rho)*(0.0  -dPhi_dn);
                }
            }
            else
            {
                if (bdry_upper_side)
                {
                    g -= a*(dt/d_rho)*(dP_dn-dPhi_dn);
                }
                else
                {
                    g += a*(dt/d_rho)*(dP_dn-dPhi_dn);
                }
            }
        }
        else
        {
            // Correct boundary conditions for the tangential velocity.
            SAMRAI::pdat::CellIndex<NDIM> i_c_intr_left = i_c_intr;
            SAMRAI::pdat::CellIndex<NDIM> i_c_bdry_left = i_c_bdry;
            i_c_intr_left(d_velocity_depth) -= 1;
            i_c_bdry_left(d_velocity_depth) -= 1;

            SAMRAI::pdat::CellIndex<NDIM> i_c_intr_rght = i_c_intr;
            SAMRAI::pdat::CellIndex<NDIM> i_c_bdry_rght = i_c_bdry;
            i_c_intr_rght(d_velocity_depth) += 1;
            i_c_bdry_rght(d_velocity_depth) += 1;

            const double   dP_dt = 0.5*(((*  P_data)(i_c_bdry_rght)+(*  P_data)(i_c_intr_rght))-((*  P_data)(i_c_bdry_left)-(*  P_data)(i_c_intr_left)))/(2.0*dx[d_velocity_depth]);
            const double dPhi_dt = 0.5*(((*Phi_data)(i_c_bdry_rght)+(*Phi_data)(i_c_intr_rght))-((*Phi_data)(i_c_bdry_left)-(*Phi_data)(i_c_intr_left)))/(2.0*dx[d_velocity_depth]);

            if (d_using_pressure_increment_form)
            {
                g -= a*(dt/d_rho)*(0.0  -dPhi_dt);
            }
            else
            {
                g -= a*(dt/d_rho)*(dP_dt-dPhi_dt);
            }
        }
    }
    return;
}// computeCorrectedBcCoefs

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
