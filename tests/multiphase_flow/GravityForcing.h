// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <SAMRAI_config.h>

// Headers for basic SAMRAI objects
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibtk/CartGridFunctionSet.h>

#ifndef included_IBAMR_multiphase_flow_GravityForcing
#define included_IBAMR_multiphase_flow_GravityForcing

/*!
 * \brief Class GravityForcing provides forcing for the momentum equations
 * due to gravity in the variable density Navier-Stokes equations
 */
class GravityForcing : public IBTK::CartGridFunction
{
public:
    GravityForcing(const std::string& object_name,
                   SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                   std::vector<double> grav_const)
        : d_object_name(object_name), d_ins_hierarchy_integrator(ins_hierarchy_integrator), d_grav_const(grav_const)
    {
        return;
    }

    bool isTimeDependent() const
    {
        return true;
    } // isTimeDependent

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > /*var*/,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 const double /*data_time*/,
                                 const bool /*initial_time*/ = false,
                                 const int coarsest_ln_in = -1,
                                 const int finest_ln_in = -1)
    {
        const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
        const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);

        // Get interpolated density variable
        const int rho_ins_idx = d_ins_hierarchy_integrator->getLinearOperatorRhoPatchDataIndex();

#if !defined(NDEBUG)
        TBOX_ASSERT(rho_ins_idx >= 0);
#endif
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& box = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > f_data = patch->getPatchData(data_idx);
                const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > rho_data =
                    patch->getPatchData(rho_ins_idx);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SAMRAI::hier::Box<NDIM>::Iterator it(SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(box, axis));
                         it;
                         it++)
                    {
                        SAMRAI::pdat::SideIndex<NDIM> s_i(it(), axis, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                        (*f_data)(s_i) = ((*rho_data)(s_i)) * d_grav_const[axis];
                    }
                }
            }
        }
        return;
    }

    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > /*var*/,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        const double /*data_time*/,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > /*patch_level*/ = nullptr)
    {
        if (initial_time)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > f_data = patch->getPatchData(data_idx);
            f_data->fillAll(0.0);
        }
    }

private:
    GravityForcing(const GravityForcing& from) = delete;
    GravityForcing& operator=(const GravityForcing& that) = delete;

    std::string d_object_name;
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_ins_hierarchy_integrator;
    std::vector<double> d_grav_const;
};

#endif
