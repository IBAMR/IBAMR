// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_vc_ins_utilities
#define included_IBAMR_vc_ins_utilities

#include <ibtk/config.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include "ibtk/CartGridFunction.h"

#include "tbox/Pointer.h"

namespace IBTK
{
class HierarchyMathOps;
}

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////
namespace IBAMR
{

namespace VcINSUtilities
{
class GravityForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    GravityForcing(const std::string& object_name,
                   SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                   std::vector<double> grav_const,
                   std::string grav_type = "FULL",
                   SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hierarchy_integrator = nullptr,
                   SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var = nullptr,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr);

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
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

    //\}

private:
    GravityForcing();

    GravityForcing(const GravityForcing& from);

    GravityForcing& operator=(const GravityForcing& that);

    std::string d_object_name;
    SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> d_ins_hierarchy_integrator = nullptr;
    std::vector<double> d_grav_const;
    std::string d_grav_type = "FULL"; // Valid options are: "FULL" and "FLOW".
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hierarchy_integrator = nullptr;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_gas_var = nullptr;
    double d_rho_neg, d_rho_pos;
    int d_num_solid_interface_cells, d_num_gas_interface_cells;
};
} // namespace VcINSUtilities

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_vc_ins_utilities
