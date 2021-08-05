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

#ifndef included_GravityForcing
#define included_GravityForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class GravityForcing provides forcing for the momentum equations
 * due to gravity in the variable density Navier-Stokes equations
 */
class GravityForcing : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    GravityForcing(const std::string& object_name,
                   Pointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                   std::vector<double> grav_const);

    /*!
     * \brief Empty destructor.
     */
    ~GravityForcing();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete GravityForcing object is
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
    GravityForcing();

    GravityForcing(const GravityForcing& from);

    GravityForcing& operator=(const GravityForcing& that);

    std::string d_object_name;
    Pointer<INSVCStaggeredHierarchyIntegrator> d_ins_hierarchy_integrator;
    std::vector<double> d_grav_const;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_GravityForcing