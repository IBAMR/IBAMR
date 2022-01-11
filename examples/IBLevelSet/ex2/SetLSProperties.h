// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_SetLSProperties
#define included_SetLSProperties

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibamr/LSInitStrategy.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * Pre processing call back functionS to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param ls_idx a patch data index for the current level set variable maintained by the integrator.
 * \param ctx is the pointer to SetLSProperties class object.
 *
 * \TODO: Let's move this out of the global namespace and use "snake case" for static function names.
 */
void callSetGasLSCallbackFunction(int ls_gas_idx,
                                  SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                  const int integrator_step,
                                  const double current_time,
                                  const bool initial_time,
                                  const bool regrid_time,
                                  void* ctx);
/*!
 * \brief Class SetLSProperties is a utility class which sets level set values on the patch hierarchy
 */
class SetLSProperties
{
public:
    /*!
     * The only constructor of this class.
     */
    SetLSProperties(const std::string& object_name, SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_gas_ops);

    /*!
     * Destructor for this class.
     */
    ~SetLSProperties() = default;

    /*!
     * Set the density based on the current solid level set information
     */
    void setLSGasPatchData(int ls_gas_idx,
                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                           const int integrator_step,
                           const double current_time,
                           const bool initial_time,
                           const bool regrid_time);

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetLSProperties() = delete;

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetLSProperties& operator=(const SetLSProperties& that) = delete;

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetLSProperties(const SetLSProperties& from) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> d_ls_gas_ops;

}; // SetLSProperties

#endif // #ifndef included_SetLSProperties
