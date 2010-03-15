// Filename: GravitationalBodyForce.C
// Last modified: <15.Mar.2010 00:22:09 griffith@griffith-macbook-pro.local>
// Created on 03 May 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "GravitationalBodyForce.h"

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
#include <CellData.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

GravitationalBodyForce::GravitationalBodyForce(
    const string& object_name,
    tbox::Pointer<tbox::Database> input_db)
    : CartGridFunction(object_name),
      d_object_name(object_name),
      d_gravitational_force(NDIM,0.0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
#endif

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    return;
}// GravitationalBodyForce

GravitationalBodyForce::~GravitationalBodyForce()
{
    // intentionally blank
    return;
}// ~GravitationalBodyForce

void
GravitationalBodyForce::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > var,
    tbox::Pointer<hier::Patch<NDIM> > patch,
    const double data_time,
    const bool initial_time,
    tbox::Pointer<hier::PatchLevel<NDIM> > level)
{
    tbox::Pointer< pdat::CellData<NDIM,double> > f_data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!f_data.isNull());
#endif
    for (int d = 0; d < NDIM; ++d)
    {
        f_data->fill(d_gravitational_force[d],d);
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
GravitationalBodyForce::getFromInput(
    tbox::Pointer<tbox::Database> db)
{
    if (!db.isNull())
    {
        if (db->keyExists("gravitational_force"))
        {
            db->getDoubleArray("gravitational_force", &d_gravitational_force[0], NDIM);
        }
    }
    return;
}// getFromInput

//////////////////////////////////////////////////////////////////////////////
