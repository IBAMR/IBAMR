// ---------------------------------------------------------------------
//
// Copyright (c) 2010 - 2018 by the IBAMR developers
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

// C++ namespace delcarations
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
    UFunction(const string& object_name, Pointer<GridGeometry<NDIM> > grid_geom, Pointer<Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~UFunction();

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
    void setDataOnPatch(int data_idx,
                        Pointer<Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > level = Pointer<PatchLevel<NDIM> >(NULL));

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
    void getFromInput(Pointer<Database> db);

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
