// Filename: SetLSProperties.h
// Created on May 24, 2018 by Nishant Nangia
//
// Copyright (c) 2002-2018, Nishant Nangia and Ameet Bhalla.
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
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param ls_idx a patch data index for the current level set variable maintained by the integrator.
 * \param ctx is the pointer to SetLSProperties class object.
 */

void callSetLSCallbackFunction(int ls_idx,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int integrator_step,
                               const double current_time,
                               const bool initial_time,
                               const bool regrid_time,
                               void* ctx);
class SetLSProperties
{
    /*!
     * \brief Class SetLSProperties is a utility class which sets level set values on the patch hierarchy
     */
public:
    /*!
     * The only constructor of this class.
     */
    SetLSProperties(const std::string& object_name, SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_ops);

    /*!
     * Destructor for this class.
     */
    ~SetLSProperties();

    /*!
     * Set the density based on the current level set information
     */
    void setLSPatchData(int ls_idx,
                        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                        const int integrator_step,
                        const double current_time,
                        const bool initial_time,
                        const bool regrid_time);

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetLSProperties();

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

    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> d_ls_ops;

}; // SetLSProperties

#endif // #ifndef included_SetLSProperties
