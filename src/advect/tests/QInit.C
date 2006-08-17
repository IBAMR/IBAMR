// Filename: QInitAndSoln.C
// Last modified: <17.Aug.2006 16:05:50 boyce@bigboy.nyconnect.com>
// Created on 19 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "QInit.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <Index.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifdef DEBUG_NO_INLINE
//#include "QInit.I"
//#endif

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

QInit::QInit(
    const string& object_name,
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
    tbox::Pointer<tbox::Database> input_db)
    : SetDataStrategy(object_name),
      d_X(NDIM)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
#endif
    d_object_name = object_name;
    d_grid_geom = grid_geom;
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_grid_geom.isNull());
#endif
    
    // Default initial values.
    const double* const XUpper = d_grid_geom->getXUpper();
    const double* const XLower = d_grid_geom->getXLower();
    
    for (int d = 0; d < NDIM; ++d)
    {
        d_X[d] = XLower[d] + 0.5*(XUpper[d] - XLower[d]);
    }

    d_init_type = "GAUSSIAN";
    
    d_gaussian_kappa = 0.01;
    
    d_zalesak_r = 0.15;
    d_zalesak_slot_w = 0.025;
    d_zalesak_slot_l = 0.1;
    
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    return;
}// QInit

QInit::~QInit()
{
    // intentionally blank
    return;
}// ~QInit

void
QInit::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > var,
    hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    (void) data_time;
    
    tbox::Pointer< pdat::CellData<NDIM,double> > Q_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!Q_data.isNull());
#endif
    const hier::Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    
    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();
    
    double r_squared;
    double X[NDIM];
    
    if (d_init_type == "GAUSSIAN")
    {
        Q_data->fillAll(0.0);
        for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            // NOTE: This assumes the lattice of Gaussians is being
            // advected in the unit square.
            int offset[NDIM];
            for (offset[0] = -2; offset[0] <= 2; ++(offset[0]))
            {
#if (NDIM>1)
                for (offset[1] = -2; offset[1] <= 2; ++(offset[1]))
                {
#endif
#if (NDIM>2)
                    for (offset[2] = -2; offset[2] <= 2; ++(offset[2]))
                    {
#endif
                        r_squared = 0.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X[d] = XLower[d] +
                                dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
                            r_squared += pow(
                                X[d]-(d_X[d]+static_cast<double>(offset[d])),2.0);
                        }
                        
                        (*Q_data)(i) +=
                            exp(-r_squared/(4.0*d_gaussian_kappa))/
                            pow(4.0*M_PI*d_gaussian_kappa,
                                0.5*static_cast<double>(NDIM));
#if (NDIM>2)
                    }
#endif
#if (NDIM>1)
                }
#endif
            }
        }
    }
    else if (d_init_type == "ZALESAK")
    {            
        for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            r_squared = 0.0;
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = XLower[d] +
                    dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
                r_squared += pow((X[d]-d_X[d]),2.0);
            }
            if ((sqrt(r_squared) > d_zalesak_r) ||
                ((abs(X[0] - d_X[0]) < d_zalesak_slot_w) &&
                 (X[1] - d_X[1]) < d_zalesak_slot_l))
            {
                (*Q_data)(i) = 0.0;
            }
            else
            {
                (*Q_data)(i) = 1.0;
            }
        }
    }
    else if (d_init_type == "SINUSOIDAL")
    {
        for (pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = XLower[d] +
                    dx[d]*(static_cast<double>(i(d)-patch_lower(d))+0.5);
            }
            (*Q_data)(i) = sin(X[0]);            
        }        
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeDataOnPatch()\n"
                   << "  invalid initialization type " << d_init_type << "\n");
    }
    return;
}// setDataOnPatch

/////////////////////////////// PRIVATE //////////////////////////////////////

void
QInit::getFromInput(
    tbox::Pointer<tbox::Database> db)
{
    if (!db.isNull())
    {
        if (db->keyExists("X"))
        {
            d_X = db->getDoubleArray("X");
        }
        
        d_init_type = db->getStringWithDefault("init_type",d_init_type);
        
        if (d_init_type == "GAUSSIAN")
        {
            d_gaussian_kappa = db->
                getDoubleWithDefault("kappa",d_gaussian_kappa);
        }
        else if (d_init_type == "ZALESAK")
        {
            d_zalesak_r = db->
                getDoubleWithDefault("zalesak_r",d_zalesak_r);
            d_zalesak_slot_w = db->
                getDoubleWithDefault("zalesak_slot_w",d_zalesak_slot_w);
            d_zalesak_slot_l = db->
                getDoubleWithDefault("zalesak_slot_l",d_zalesak_slot_l);
        }
        else if (d_init_type == "SINUSOIDAL")
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

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the QInit class.
///
//////////////////////////////////////////////////////////////////////

template class SAMRAI::tbox::Pointer<QInit>;

//////////////////////////////////////////////////////////////////////////////
