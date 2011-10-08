// Filename: FACPreconditioner.I
// Created on 16 Sep 2003 by Boyce Griffith
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

inline void
FACPreconditioner::setInitialGuessNonzero(
    bool /*initial_guess_nonzero*/)
{
    // intentionally blank
    return;
}// setInitialGuessNonzero

inline bool
FACPreconditioner::getInitialGuessNonzero() const
{
    // intentionally blank
    return false;
}// getInitialGuessNonzero

inline void
FACPreconditioner::setMaxIterations(
    int max_iterations)
{
    if (max_iterations > 1)
    {
        TBOX_ERROR(d_object_name << "::setMaxIterations()\n"
                   << "  class IBTK::FACPreconditioner only performs a single iteration" << std::endl);
    }
    return;
}// setMaxIterations

inline int
FACPreconditioner::getMaxIterations() const
{
    // intentionally blank
    return 1;
}// getMaxIterations

inline void
FACPreconditioner::setAbsoluteTolerance(
    double /*abs_residual_tol*/)
{
    // intentionally blank
    return;
}//setAbsoluteTolerance

inline double
FACPreconditioner::getAbsoluteTolerance() const
{
    // intentionally blank
    return 0.0;
}// getAbsoluteTolerance

inline void
FACPreconditioner::setRelativeTolerance(
    double /*rel_residual_tol*/)
{
    // intentionally blank
    return;
}//setRelativeTolerance

inline double
FACPreconditioner::getRelativeTolerance() const
{
    // intentionally blank
    return 0.0;
}// getRelativeTolerance

inline void
FACPreconditioner::setMGCycleType(
    MGCycleType cycle_type)
{
    d_cycle_type = cycle_type;
    return;
}// setMGCycleType

inline MGCycleType
FACPreconditioner::getMGCycleType() const
{
    return d_cycle_type;
}// getMGCycleType

inline void
FACPreconditioner::setNumPreSmoothingSweeps(
    int num_pre_sweeps)
{
    d_num_pre_sweeps = num_pre_sweeps;
    return;
}// setNumPreSmoothingSweeps

inline int
FACPreconditioner::getNumPreSmoothingSweeps() const
{
    return d_num_pre_sweeps;
}// getNumPreSmoothingSweeps

inline void
FACPreconditioner::setNumPostSmoothingSweeps(
    int num_post_sweeps)
{
    d_num_post_sweeps = num_post_sweeps;
    return;
}// setNumPostSmoothingSweeps

inline int
FACPreconditioner::getNumPostSmoothingSweeps() const
{
    return d_num_post_sweeps;
}// getNumPostSmoothingSweeps

inline int
FACPreconditioner::getNumIterations() const
{
    // intentionally blank
    return 1;
}// getNumIterations

inline double
FACPreconditioner::getResidualNorm() const
{
    // intentionally blank
    return 0.0;
}// getResidualNorm

inline void
FACPreconditioner::enableLogging(
    bool enabled)
{
    d_do_log = enabled;
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
