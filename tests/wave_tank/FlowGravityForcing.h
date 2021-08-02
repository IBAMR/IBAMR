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

#ifndef included_FlowGravityForcing
#define included_FlowGravityForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class FlowGravityForcing provides forcing for the momentum equations
 * due to gravity in the variable coefficient Navier-Stokes equations, reconstructed
 * from the flow density field.
 */
class FlowGravityForcing : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    FlowGravityForcing(const std::string& object_name,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_hierarchy_integrator,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var,
                       std::vector<double> grav_const);

    /*!
     * \brief Empty destructor.
     */
    ~FlowGravityForcing();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete FlowGravityForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 Pointer<hier::Variable<NDIM> > var,
                                 Pointer<PatchHierarchy<NDIM> > hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<hier::Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL));

    //\}

private:
    FlowGravityForcing();

    FlowGravityForcing(const FlowGravityForcing& from);

    FlowGravityForcing& operator=(const FlowGravityForcing& that);

    std::string d_object_name;
    Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hierarchy_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_gas_var;
    double d_rho_neg, d_rho_pos;
    std::vector<double> d_grav_const;
    int d_num_solid_interface_cells, d_num_gas_interface_cells;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FlowGravityForcing
