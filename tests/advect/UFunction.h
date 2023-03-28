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

#ifndef included_IBAMR_advect_UFunction
#define included_IBAMR_advect_UFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Method to initialize the value of the advection velocity u.
 */
class UFunction : public CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    UFunction(const std::string& object_name,
              SAMRAI::tbox::Pointer<SAMRAI::hier::GridGeometry<NDIM> > grid_geom,
              SAMRAI::tbox::Pointer<Database> input_db)
        : CartGridFunction(object_name),
          d_object_name(object_name),
          d_grid_geom(grid_geom),
          d_X(),
          d_init_type("UNIFORM"),
          d_uniform_u()
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(!object_name.empty());
        TBOX_ASSERT(grid_geom);
#endif

        // Default initial values.
        const double* const x_upper = d_grid_geom->getXUpper();
        const double* const x_lower = d_grid_geom->getXLower();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_X[d] = x_lower[d] + 0.5 * (x_upper[d] - x_lower[d]);
            d_uniform_u[d] = 1.0;
        }

        // Initialize object with data read from the input database.
        getFromInput(input_db);

        return;
    } // UFunction

    /*!
     * Indicates whether the concrete CartGridFunction object is time dependent.
     */
    bool isTimeDependent() const
    {
        return true;
    }

    /*!
     * Set the data on the patch interior to some initial values.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > /*var*/,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        const double /*data_time*/,
                        const bool /*initial_time*/,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > /*level*/)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > u_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(u_data);
#endif

        if (d_init_type == "UNIFORM")
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                u_data->getArrayData(axis).fillAll(d_uniform_u[axis]);
            }
        }
        else if (d_init_type == "VORTEX")
        {
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();

            IBTK::VectorNd X;

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (SAMRAI::pdat::FaceIterator<NDIM> it(patch_box, axis); it; it++)
                {
                    const SAMRAI::pdat::FaceIndex<NDIM>& i = it();
                    const SAMRAI::hier::Index<NDIM>& cell_idx = i.toCell(1);

                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] +
                               dx[d] * (static_cast<double>(cell_idx(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
                    }

                    // 2D vortex
                    if (axis == 0)
                    {
                        (*u_data)(i) = (X[1] - d_X[axis]);
                    }
                    else if (axis == 1)
                    {
                        (*u_data)(i) = (d_X[axis] - X[0]);
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
            TBOX_ERROR(d_object_name << "::setDataOnPatch()\n"
                                     << "  invalid initialization type " << d_init_type << "\n");
        }
        return;
    } // setDataOnPatch

protected:
private:
    UFunction(const UFunction& from) = delete;
    UFunction& operator=(const UFunction& that) = delete;

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

            if (d_init_type == "UNIFORM")
            {
                if (db->keyExists("uniform_u"))
                {
                    db->getDoubleArray("uniform_u", d_uniform_u.data(), NDIM);
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
    } // getFromInput

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    std::string d_object_name;

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
    std::string d_init_type;

    /*
     * Parameters for uniform constant velocity.
     */
    IBTK::VectorNd d_uniform_u;
};

#endif //#ifndef included_UFunction
