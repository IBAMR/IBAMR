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

#ifndef included_VelocityInitialCondition
#define included_VelocityInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>

// Application includes
#include "LSLocateFluidInterface.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class VelocityInitialCondition provides forcing an initial condition for velocity
 * based on the initial level set information.
 */
class VelocityInitialCondition : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    VelocityInitialCondition(const std::string& object_name,
                             const double num_interface_cells,
                             std::vector<double> inside_velocity,
                             std::vector<double> outside_velocity,
                             CircularInterface init_circle);

    /*!
     * \brief Empty destructor.
     */
    ~VelocityInitialCondition();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete VelocityInitialCondition object is
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
    VelocityInitialCondition();

    VelocityInitialCondition(const VelocityInitialCondition& from);

    VelocityInitialCondition& operator=(const VelocityInitialCondition& that);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Number of interface cells over which to smooth the material properties
     */
    double d_num_interface_cells;

    /*!
     * Velocities of the fluid level set
     */
    std::vector<double> d_inside_velocity;
    std::vector<double> d_outside_velocity;

    /*!
     * Initial level set information.
     */
    CircularInterface d_init_circle;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VelocityInitialCondition
