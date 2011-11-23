// Filename: PETScSNESJacobianJOWrapper.I
// Created on 23 Aug 2006 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
PETScSNESJacobianJOWrapper::PETScSNESJacobianJOWrapper(
    const std::string& object_name,
    const SNES& petsc_snes,
    PetscErrorCode (* const petsc_snes_form_jac)(SNES,Vec,Mat*,Mat*,MatStructure*,void*),
    void* const petsc_snes_jac_ctx)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_do_log(false),
      d_petsc_snes(petsc_snes),
      d_petsc_snes_jac(PETSC_NULL),
      d_petsc_snes_form_jac(petsc_snes_form_jac),
      d_petsc_snes_jac_ctx(petsc_snes_jac_ctx),
      d_x(NULL),
      d_y(NULL),
      d_z(NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_y(PETSC_NULL),
      d_petsc_z(PETSC_NULL)
{
    // intentionally blank
    return;
}// PETScSNESJacobianJOWrapper()

inline
PETScSNESJacobianJOWrapper::~PETScSNESJacobianJOWrapper()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~PETScSNESJacobianJOWrapper()

inline const SNES&
PETScSNESJacobianJOWrapper::getPETScSNES() const
{
    return d_petsc_snes;
}// getPETScSNES

inline PetscErrorCode
(*PETScSNESJacobianJOWrapper::getPETScSNESFormJacobian())(SNES,Vec,Mat*,Mat*,MatStructure*,void*)
{
    return d_petsc_snes_form_jac;
}// getPETScSNESFormJacobian

inline void*
PETScSNESJacobianJOWrapper::getPETScSNESJacobianContext() const
{
    return d_petsc_snes_jac_ctx;
}// getPETScSNESJacobianContext

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
