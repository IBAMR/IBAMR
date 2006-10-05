// Filename: USet.C
// Last modified: <04.Sep.2006 02:15:22 boyce@bigboy.nyconnect.com>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "USet.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

#include <ArrayData.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <FaceData.h>
#include <FaceIndex.h>
#include <FaceIterator.h>
#include <Index.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

USet::USet(
    const string& object_name,
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
    tbox::Pointer<tbox::Database> input_db)
    : SetDataStrategy(object_name),
      d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_X(NDIM),
      d_init_type("UNIFORM"),
      d_uniform_u(NDIM)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
#endif

    // Default initial values.
    const double* const XUpper = d_grid_geom->getXUpper();
    const double* const XLower = d_grid_geom->getXLower();
    for (int d = 0; d < NDIM; ++d)
    {
        d_X[d] = XLower[d] + 0.5*(XUpper[d] - XLower[d]);
        d_uniform_u[d] = 1.0;
    }

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    return;
}// USet

USet::~USet()
{
    // intentionally blank
    return;
}// ~USet

void
USet::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > var,
    hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    if (initial_time)
    {
        tbox::Pointer< pdat::FaceData<NDIM,double> > u_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!u_data.isNull());
#endif

        if (d_init_type == "UNIFORM")
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                u_data->getArrayData(axis).fillAll(d_uniform_u[axis]);
            }
        }
        else if (d_init_type == "VORTEX")
        {
            const hier::Box<NDIM>& patch_box = patch.getBox();
            const hier::Index<NDIM>& patch_lower = patch_box.lower();
            tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

            const double* const XLower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();

            double X[NDIM];

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (pdat::FaceIterator<NDIM> it(patch_box,axis); it; it++)
                {
                    const pdat::FaceIndex<NDIM>& i = it();
                    const hier::Index<NDIM>& cell_idx = i.toCell(1);

                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (d != axis)
                        {
                            X[d] =
                                XLower[d] +
                                dx[d]*(static_cast<double>(cell_idx(d)-patch_lower(d))+0.5);
                        }
                        else
                        {
                            X[d] =
                                XLower[d] +
                                dx[d]*(static_cast<double>(cell_idx(d)-patch_lower(d)));
                        }
                    }

                    // 2D vortex
                    if (axis == 0)
                    {
                        (*u_data)(i) = X[1] - d_X[axis];
                    }
                    else if (axis == 1)
                    {
                        (*u_data)(i) = d_X[axis] - X[0];
                    }
                    else
                    {
                        (*u_data)(i) = 0.0;
                    }
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeDataOnPatch()\n"
                       << "  invalid initialization type " << d_init_type << "\n");
        }
    }
    return;
}// setDataOnPatch

/////////////////////////////// PRIVATE //////////////////////////////////////

void
USet::getFromInput(
    tbox::Pointer<tbox::Database> db)
{
    if (!db.isNull())
    {
        if (db->keyExists("X"))
        {
            d_X = db->getDoubleArray("X");
        }

        d_init_type = db->getStringWithDefault("init_type",d_init_type);

        if (d_init_type == "UNIFORM")
        {
            if (db->keyExists("uniform_u"))
            {
                d_uniform_u = db->getDoubleArray("uniform_u");
            }
        }
        else if (d_init_type == "VORTEX")
        {
            // intentionally blank
        }
        else
        {
            TBOX_ERROR(d_object_name << "::getFromInput()\n"
                       << "  invalid initialization type " << d_init_type << "\n");
        }
    }
    return;
}// getFromInput

//////////////////////////////////////////////////////////////////////////////
