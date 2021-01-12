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
#include <ibamr/LSInitStrategy.h>

#ifndef included_IBAMR_multiphase_flow_SetLSProperties
#define included_IBAMR_multiphase_flow_SetLSProperties

class SetLSProperties
{
    /*!
     * \brief Class SetLSProperties is a utility class which sets level set values on the patch hierarchy
     */
public:
    /*!
     * The only constructor of this class.
     */
    SetLSProperties(const std::string& object_name,
                    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_solid_ops,
                    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_gas_ops)
        : d_object_name(object_name), d_ls_solid_ops(ls_solid_ops), d_ls_gas_ops(ls_gas_ops)
    {
        // intentionally left blank
        return;
    } // SetLSProperties

    /*!
     * Set the density based on the current solid level set information
     */
    void setLSSolidPatchData(int ls_solid_idx,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int integrator_step,
                             const double current_time,
                             const bool initial_time,
                             const bool regrid_time)
    {
        // If at the regrid time, force reinitialization
        TBOX_ASSERT(d_ls_solid_ops);
        d_ls_solid_ops->setReinitializeLSData(regrid_time);
        d_ls_solid_ops->initializeLSData(ls_solid_idx, hier_math_ops, integrator_step, current_time, initial_time);

        return;
    } // setLSSolidPatchData

    /*!
     * Set the density based on the current solid level set information
     */
    void setLSGasPatchData(int ls_gas_idx,
                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                           const int integrator_step,
                           const double current_time,
                           const bool initial_time,
                           const bool regrid_time)
    {
        // If at the regrid time, force reinitialization
        TBOX_ASSERT(d_ls_gas_ops);
        d_ls_gas_ops->setReinitializeLSData(regrid_time);
        d_ls_gas_ops->initializeLSData(ls_gas_idx, hier_math_ops, integrator_step, current_time, initial_time);

        return;
    } // setLSSolidPatchData

private:
    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetLSProperties& operator=(const SetLSProperties& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetLSProperties(const SetLSProperties& from);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> d_ls_solid_ops;
    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> d_ls_gas_ops;

}; // SetLSProperties

inline void
callSetSolidLSCallbackFunction(int ls_solid_idx,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int integrator_step,
                               const double current_time,
                               const bool initial_time,
                               const bool regrid_time,
                               void* ctx)
{
    // Set the density from the level set information
    static SetLSProperties* ptr_SetLSProperties = static_cast<SetLSProperties*>(ctx);
    ptr_SetLSProperties->setLSSolidPatchData(
        ls_solid_idx, hier_math_ops, integrator_step, current_time, initial_time, regrid_time);

    return;
} // callSetSolidLSCallbackFunction

inline void
callSetGasLSCallbackFunction(int ls_gas_idx,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int integrator_step,
                             const double current_time,
                             const bool initial_time,
                             const bool regrid_time,
                             void* ctx)
{
    // Set the density from the level set information
    static SetLSProperties* ptr_SetLSProperties = static_cast<SetLSProperties*>(ctx);
    ptr_SetLSProperties->setLSGasPatchData(
        ls_gas_idx, hier_math_ops, integrator_step, current_time, initial_time, regrid_time);

    return;

} // callSetGasLSCallbackFunction

#endif
