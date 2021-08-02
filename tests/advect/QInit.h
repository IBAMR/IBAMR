// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_QInit
#define included_QInit

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

#include <ibtk/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the advected scalar Q.
 */
class QInit : public CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    QInit(const string& object_name,
          SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom,
          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
        : CartGridFunction(object_name),
          d_object_name(object_name),
          d_grid_geom(grid_geom),
          d_X(),
          d_init_type("GAUSSIAN"),
          d_gaussian_kappa(0.01),
          d_zalesak_r(0.15),
          d_zalesak_slot_w(0.025),
          d_zalesak_slot_l(0.1)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(!object_name.empty());
        TBOX_ASSERT(grid_geom);
#endif
        d_object_name = object_name;
        d_grid_geom = grid_geom;
#if !defined(NDEBUG)
        TBOX_ASSERT(d_grid_geom);
#endif

        // Default initial values.
        const double* const x_upper = d_grid_geom->getXUpper();
        const double* const x_lower = d_grid_geom->getXLower();

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_X[d] = x_lower[d] + 0.5 * (x_upper[d] - x_lower[d]);
        }

        d_init_type = "GAUSSIAN";

        d_gaussian_kappa = 0.01;

        d_zalesak_r = 0.15;
        d_zalesak_slot_w = 0.025;
        d_zalesak_slot_l = 0.1;

        // Initialize object with data read from the input database.
        getFromInput(input_db);

        return;
    } // QInit

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool isTimeDependent() const
    {
        return true;
    }

    /*!
     * Set the data on the patch interior to the exact answer.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<Variable<NDIM> > /*var*/,
                        SAMRAI::tbox::Pointer<Patch<NDIM> > patch,
                        const double data_time,
                        const bool /*initial_time*/,
                        SAMRAI::tbox::Pointer<PatchLevel<NDIM> > /*level*/)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > Q_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(Q_data);
#endif
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

        const double* const x_lower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();

        double r_squared;
        IBTK::VectorNd X;
        const double t = data_time;

        Q_data->fillAll(0.0);

        if (d_init_type == "GAUSSIAN")
        {
            for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const SAMRAI::hier::Index<NDIM>& i = ic();
                // NOTE: This assumes the lattice of Gaussians are being advected
                // and diffused in the unit square.
                std::array<int, NDIM> offset;
                for (offset[0] = -2; offset[0] <= 2; ++(offset[0]))
                {
                    for (offset[1] = -2; offset[1] <= 2; ++(offset[1]))
                    {
#if (NDIM > 2)
                        for (offset[2] = -2; offset[2] <= 2; ++(offset[2]))
                        {
#endif
                            r_squared = 0.0;
                            for (unsigned int d = 0; d < NDIM; ++d)
                            {
                                X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                                r_squared += pow(X[d] - (d_X[d] + static_cast<double>(offset[d])), 2.0);
                            }
                            (*Q_data)(i) +=
                                exp(-r_squared / (4.0 * d_gaussian_kappa)) /
                                pow(4.0 * M_PI * d_gaussian_kappa * (1.0 + t), 0.5 * static_cast<double>(NDIM));
#if (NDIM > 2)
                        }
#endif
                    }
                }
            }
        }
        else if (d_init_type == "ZALESAK")
        {
            for (SAMRAI::pdat::CellIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const SAMRAI::hier::Index<NDIM>& i = ic();
                r_squared = 0.0;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                    r_squared += pow((X[d] - d_X[d]), 2.0);
                }
                if ((sqrt(r_squared) > d_zalesak_r) ||
                    ((abs(X[0] - d_X[0]) < d_zalesak_slot_w) && (X[1] - d_X[1]) < d_zalesak_slot_l))
                {
                    (*Q_data)(i) = 0.0;
                }
                else
                {
                    (*Q_data)(i) = 1.0;
                }
            }
        }
        else
        {
            TBOX_ERROR(d_object_name << "::initializeDataOnPatch()\n"
                                     << "  invalid initialization type " << d_init_type << "\n");
        }
        return;
    } // setDataOnPatch

protected:
private:
    QInit(const QInit& from) = delete;
    QInit& operator=(const QInit& that) = delete;

    /*!
     * Read input values, indicated above, from given database.
     */
    void getFromInput(Pointer<Database> db)
    {
        if (db)
        {
            if (db->keyExists("X"))
            {
                db->getDoubleArray("X", d_X.data(), NDIM);
            }

            d_init_type = db->getStringWithDefault("init_type", d_init_type);

            if (d_init_type == "GAUSSIAN")
            {
                d_gaussian_kappa = db->getDoubleWithDefault("kappa", d_gaussian_kappa);
            }
            else if (d_init_type == "ZALESAK")
            {
                d_zalesak_r = db->getDoubleWithDefault("zalesak_r", d_zalesak_r);
                d_zalesak_slot_w = db->getDoubleWithDefault("zalesak_slot_w", d_zalesak_slot_w);
                d_zalesak_slot_l = db->getDoubleWithDefault("zalesak_slot_l", d_zalesak_slot_l);
            }
            else
            {
                TBOX_ERROR(d_object_name << "::getFromInput()\n"
                                         << "  invalid initialization type " << d_init_type << "\n");
            }
        }
        return;
    } // getFromInput

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    string d_object_name;

    /*
     * The grid geometry.
     */
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*
     * The center of the initial data.
     */
    IBTK::VectorNd d_X;

    /*
     * The initialization type.
     */
    string d_init_type;

    /*
     * Parameters for Gaussian initial conditions.
     */
    double d_gaussian_kappa;

    /*
     * Parameters for the Zalesak slotted cylinder.
     */
    double d_zalesak_r;
    double d_zalesak_slot_w;
    double d_zalesak_slot_l;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "QInit.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_QInit
