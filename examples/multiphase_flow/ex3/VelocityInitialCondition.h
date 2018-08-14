// Filename: VelocityInitialCondition.h
// Created on 27 June 2018 by Nishant Nangia
//
// Copyright (c) 2002-2018, Nishant Nangia
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_VelocityInitialCondition
#define included_VelocityInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>

// Application includes
#include "LSLocateCircularInterface.h"

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
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 Pointer<Variable<NDIM> > var,
                                 Pointer<PatchHierarchy<NDIM> > hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1);

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
