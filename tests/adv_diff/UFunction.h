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

#ifndef included_UFunction
#define included_UFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

// C++ namespace declarations.
#include <ibamr/app_namespaces.h>

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
    UFunction(const string& object_name, Pointer<GridGeometry<NDIM> > grid_geom, Pointer<Database> input_db)
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
     * \brief Destructor.
     */
    ~UFunction() = default;

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
                        Pointer<Variable<NDIM> > /*var*/,
                        Pointer<Patch<NDIM> > patch,
                        const double /*data_time*/,
                        const bool /*initial_time*/,
                        Pointer<PatchLevel<NDIM> > /*level*/ = Pointer<PatchLevel<NDIM> >(NULL))
    {
        Pointer<FaceData<NDIM, double> > u_data = patch->getPatchData(data_idx);
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
            const Box<NDIM>& patch_box = patch->getBox();
            const hier::Index<NDIM>& patch_lower = patch_box.lower();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();

            VectorNd X;

            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (FaceIterator<NDIM> it(patch_box, axis); it; it++)
                {
                    const FaceIndex<NDIM>& i = it();
                    const hier::Index<NDIM>& cell_idx = i.toCell(1);

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
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    UFunction();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    UFunction(const UFunction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    UFunction& operator=(const UFunction& that);

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
    string d_object_name;

    /*
     * The grid geometry.
     */
    Pointer<CartesianGridGeometry<NDIM> > d_grid_geom;

    /*
     * The center of the initial data.
     */
    VectorNd d_X;

    /*
     * The initialization type.
     */
    string d_init_type;

    /*
     * Parameters for uniform constant velocity.
     */
    VectorNd d_uniform_u;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "UFunction.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_UFunction
