// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_LevelSetSolidInitialCondition
#define included_LevelSetSolidInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>

// Application includes
#include "LSLocateTrapezoidalInterface.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class LevelSetSolidInitialCondition provides an initial condition for
 * the level set function.
 */
class LevelSetSolidInitialCondition : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LevelSetSolidInitialCondition(const std::string& object_name, TrapezoidalInterface* trapezoid);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetSolidInitialCondition();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete LevelSetSolidInitialCondition object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL));

    //\}

private:
    LevelSetSolidInitialCondition();

    LevelSetSolidInitialCondition(const LevelSetSolidInitialCondition& from);

    LevelSetSolidInitialCondition& operator=(const LevelSetSolidInitialCondition& that);

    /*!
     * Point to line distance
     */
    double pointToLineDistance(IBTK::Vector X0, IBTK::Vector P1, IBTK::Vector P2);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Initial level set information.
     */
    TrapezoidalInterface* d_init_trapezoid;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LevelSetSolidInitialCondition
