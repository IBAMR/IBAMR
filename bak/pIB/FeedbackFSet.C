// Filename: FeedbackFSet.C
// Last modified: <01.Feb.2007 23:41:47 boyce@bigboy.nyconnect.com>
// Created on 20 Nov 2006 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "FeedbackFSet.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <ArrayData.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <Index.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FeedbackFSet::FeedbackFSet(
    const string& object_name,
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
    tbox::Pointer<tbox::Database> input_db)
    : SetDataStrategy(object_name),
      d_U_var(NULL),
      d_U_context(NULL),
      d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_inflow_vel(-300.0),
      d_alpha(1.0e3),
      d_tau(5.0e-4)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
#endif
    if (!input_db.isNull())
    {
        d_inflow_vel = input_db->getDoubleWithDefault("inflow_vel",d_inflow_vel);
        d_alpha = input_db->getDoubleWithDefault("alpha",d_alpha);
        d_tau = input_db->getDoubleWithDefault("tau",d_tau);
    }
    return;
}// FeedbackFSet

FeedbackFSet::~FeedbackFSet()
{
    // intentionally blank
    return;
}// ~FeedbackFSet

void
FeedbackFSet::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > var,
    hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    tbox::Pointer<pdat::CellData<NDIM,double> > F_data = patch.getPatchData(data_idx);
    tbox::Pointer<pdat::CellData<NDIM,double> > U_data = patch.getPatchData(d_U_var, d_U_context);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!F_data.isNull());
    assert(!U_data.isNull());
#endif

    // Supply a feedback force to approximately impose the desired
    // velocities at the periodic boundary.
    if (!d_grid_geom->getDomainIsSingleBox()) TBOX_ERROR("physical domain must be a single box...\n");

    const hier::Box<NDIM>& domain_box = d_grid_geom->getPhysicalDomain()(0);

    hier::Box<NDIM> lower_box = domain_box;
    hier::Box<NDIM> upper_box = domain_box;
    lower_box.upper(1) = domain_box.lower(1);
    upper_box.lower(1) = domain_box.upper(1);

    const hier::IntVector<NDIM>& ratio = patch.getPatchGeometry()->getRatio();

    const hier::Box<NDIM> refined_lower_box =
        hier::Box<NDIM>::refine(lower_box,ratio);
    const hier::Box<NDIM> refined_upper_box =
        hier::Box<NDIM>::refine(upper_box,ratio);

    const hier::Box<NDIM>& patch_box = patch.getBox();

    const double t_fac = 1.0 - exp(-data_time/d_tau);

    F_data->fillAll(0.0);

    for (pdat::CellIterator<NDIM> ic(patch_box*refined_lower_box); ic; ic++)
    {
        const hier::Index<NDIM>& i = ic();
        (*F_data)(i,0) = d_alpha*(0.0               -(*U_data)(i,0));
        (*F_data)(i,1) = d_alpha*(t_fac*d_inflow_vel-(*U_data)(i,1));
    }

    for (pdat::CellIterator<NDIM> ic(patch_box*refined_upper_box); ic; ic++)
    {
        const hier::Index<NDIM>& i = ic();
        (*F_data)(i,0) = d_alpha*(0.0               -(*U_data)(i,0));
        (*F_data)(i,1) = d_alpha*(t_fac*d_inflow_vel-(*U_data)(i,1));
    }

    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
