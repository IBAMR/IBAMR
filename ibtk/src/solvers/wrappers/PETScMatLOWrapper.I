// Filename: PETScMatLOWrapper.I
// Created on 26 Dec 2003 by Boyce Griffith
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

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
PETScMatLOWrapper::PETScMatLOWrapper(
    const std::string& object_name,
    const Mat& petsc_mat)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_do_log(false),
      d_petsc_mat(petsc_mat),
      d_x(NULL),
      d_y(NULL),
      d_z(NULL),
      d_petsc_x(PETSC_NULL),
      d_petsc_y(PETSC_NULL),
      d_petsc_z(PETSC_NULL)
{
    // intentionally blank
    return;
}// PETScMatLOWrapper()

inline
PETScMatLOWrapper::~PETScMatLOWrapper()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~PETScMatLOWrapper()

inline bool
PETScMatLOWrapper::isSymmetric() const
{
    PetscTruth set, is_symmetric;
    int ierr = MatIsSymmetricKnown(d_petsc_mat,&set,&is_symmetric); IBTK_CHKERRQ(ierr);
    return (set && is_symmetric);
}// isSymmetric

inline const Mat&
PETScMatLOWrapper::getPETScMat() const
{
    return d_petsc_mat;
}// getPETScMat

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
