// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_LevelSetInitialConditionHexagram
#define included_LevelSetInitialConditionHexagram

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAIVariable.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class LevelSetInitialCondition provides an initial condition for
 * the level set function.
 */
class LevelSetInitialConditionHexagram : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LevelSetInitialConditionHexagram(const std::string& object_name,
                                     const SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geom,
                                     const IBTK::VectorNd& origin);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetInitialConditionHexagram() = default;

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete LevelSetInitialCondition object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Indicates whether the LevelSetInitialCondition corresponds to inner
     * or outer cylinder.
     */
    void setCylinderType(std::string type);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<SAMRAIVariable> var,
                        Pointer<SAMRAIPatch> patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<SAMRAIPatchLevel> patch_level = nullptr) override;

    //\}

private:
    /*!
     * Deleted default constructor.
     */
    LevelSetInitialConditionHexagram() = delete;

    /*!
     * Deleted copy constructor.
     */
    LevelSetInitialConditionHexagram(const LevelSetInitialConditionHexagram& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    LevelSetInitialConditionHexagram& operator=(const LevelSetInitialConditionHexagram& that) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * The Cartesian grid geometry object provides the extents of the
     * computational domain.
     */
    SAMRAIPointer<SAMRAICartesianGridGeometry> d_grid_geom;

    /*!
     * Origin of the geometry.
     */
    IBTK::VectorNd d_origin;
};

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_LevelSetInitialCondition
