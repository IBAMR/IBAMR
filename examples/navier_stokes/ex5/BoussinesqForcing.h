// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_BoussinesqForcing
#define included_BoussinesqForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIVariable.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class BoussinesqForcing provides forcing for the momentum equations
 * based on the Boussinesq approximation to the variable-density incompressible
 * Navier-Stokes equations.
 */
class BoussinesqForcing : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    BoussinesqForcing(Pointer<SAMRAIVariable> T_var,
                      Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
                      int gamma);

    /*!
     * \brief Empty destructor.
     */
    ~BoussinesqForcing();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete BoussinesqForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 Pointer<SAMRAIVariable> var,
                                 Pointer<SAMRAIPatchHierarchy> hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<SAMRAIVariable> var,
                        Pointer<SAMRAIPatch> patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<SAMRAIPatchLevel> patch_level = nullptr);

    //\}

private:
    BoussinesqForcing();

    BoussinesqForcing(const BoussinesqForcing& from);

    BoussinesqForcing& operator=(const BoussinesqForcing& that);

    Pointer<SAMRAIVariable> d_T_var;
    Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;
    double d_gamma;
};

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_BoussinesqForcing
