// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2018 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_LevelSetInitialCondition
#define included_LevelSetInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/app_namespaces.h>

#include <ibtk/muParserCartGridFunction.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class LevelSetInitialCondition provides an initial condition for
 * the level set function.
 */
class LevelSetInitialCondition : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LevelSetInitialCondition(const std::string& object_name,
                             const double radius,
                             const IBTK::VectorNd& origin,
                             const bool fluid_is_interior_to_cylinder);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetInitialCondition();

    /*!
     * \brief Indicates whether the concrete LevelSetInitialCondition object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL));

    //\}

private:
    LevelSetInitialCondition();

    LevelSetInitialCondition(const LevelSetInitialCondition& from);

    LevelSetInitialCondition& operator=(const LevelSetInitialCondition& that);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Radius of cylinder.
     */
    double d_radius;

    /*!
     * Origin of cylinder.
     */
    IBTK::VectorNd d_origin;

    /*!
     * Boolean to identify whether fluid is considered inside the cylinder or outside.
     * Default is set false to denote the annulus case.
     */
    bool d_fluid_is_interior_to_cylinder = false;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LevelSetInitialCondition
