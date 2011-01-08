// Filename: PETScPCLSWrapper.I
// Created on 19 Oct 2003 by Boyce Griffith
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

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
PETScPCLSWrapper::PETScPCLSWrapper(
    const std::string& object_name,
    const PC& petsc_pc)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_do_log(false),
      d_petsc_pc(petsc_pc),
      d_x(NULL),
      d_b(NULL),
      d_petsc_x(static_cast<Vec>(NULL)),
      d_petsc_b(static_cast<Vec>(NULL))
{
    // intentionally blank
    return;
}// PETScPCLSWrapper()

inline
PETScPCLSWrapper::~PETScPCLSWrapper()
{
    if (d_is_initialized) deallocateSolverState();
    return;
}// ~PETScPCLSWrapper()

inline const PC&
PETScPCLSWrapper::getPETScPC() const
{
    return d_petsc_pc;
}// getPETScPC

inline void
PETScPCLSWrapper::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    TBOX_WARNING("PETScPCLSWrapper::setInitialGuessNonzero() not supported" << std::endl);
    return;
}// setInitialGuessNonzero

inline bool
PETScPCLSWrapper::getInitialGuessNonzero() const
{
    TBOX_WARNING("PETScPCLSWrapper::getInitialGuessNonzero() not supported" << std::endl);
    return true;
}// getInitialGuessNonzero

inline void
PETScPCLSWrapper::setMaxIterations(
    int max_iterations)
{
    TBOX_WARNING("PETScPCLSWrapper::setMaxIterations() not supported" << std::endl);
    return;
}// setMaxIterations

inline int
PETScPCLSWrapper::getMaxIterations() const
{
    TBOX_WARNING("PETScPCLSWrapper::getMaxIterations() not supported" << std::endl);
    return 0;
}// getMaxIterations

inline void
PETScPCLSWrapper::setAbsoluteTolerance(
    double abs_residual_tol)
{
    TBOX_WARNING("PETScPCLSWrapper::setAbsoluteTolerance() not supported" << std::endl);
    return;
}//setAbsoluteTolerance

inline double
PETScPCLSWrapper::getAbsoluteTolerance() const
{
    TBOX_WARNING("PETScPCLSWrapper::getAbsoluteTolerance() not supported" << std::endl);
    return 0.0;
}// getAbsoluteTolerance

inline void
PETScPCLSWrapper::setRelativeTolerance(
    double rel_residual_tol)
{
    TBOX_WARNING("PETScPCLSWrapper::setRelativeTolerance() not supported" << std::endl);
    return;
}//setRelativeTolerance

inline double
PETScPCLSWrapper::getRelativeTolerance() const
{
    TBOX_WARNING("PETScPCLSWrapper::getRelativeTolerance() not supported" << std::endl);
    return 0.0;
}// getRelativeTolerance

inline int
PETScPCLSWrapper::getNumIterations() const
{
    TBOX_WARNING("PETScPCLSWrapper::getNumIteratons() not supported" << std::endl);
    return 0;
}// getNumIterations

inline double
PETScPCLSWrapper::getResidualNorm() const
{
    TBOX_WARNING("PETScPCLSWrapper::getResidualNorm() not supported" << std::endl);
    return 0.0;
}// getResidualNorm

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
