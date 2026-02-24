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

#ifndef included_LevelSetInitialCondition
#define included_LevelSetInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIPatch.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIVariable.h"

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

// Struct to maintain the properties of the rectangular interface
struct RectangularInterface
{
    Eigen::Vector3d X0;
    double W, H;
    double rho_solid;
    double g_y;
};

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
    LevelSetInitialCondition(const std::string& object_name, RectangularInterface init_rectangle);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetInitialCondition();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete LevelSetInitialCondition object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAIPointer<SAMRAIVariable> var,
                        SAMRAIPointer<SAMRAIPatch> patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAIPointer<SAMRAIPatchLevel> patch_level = SAMRAIPointer<SAMRAIPatchLevel>(nullptr));

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
     * Initial level set information.
     */
    RectangularInterface d_init_rectangle;
};

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_LevelSetInitialCondition
