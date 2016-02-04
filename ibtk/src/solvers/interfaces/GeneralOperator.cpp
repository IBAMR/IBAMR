// Filename: GeneralOperator.cpp
// Created on 18 Nov 2003 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <limits>
#include <ostream>
#include <string>
#include <utility>

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/GeneralOperator.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

GeneralOperator::GeneralOperator(const std::string& object_name, bool homogeneous_bc)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_homogeneous_bc(homogeneous_bc),
      d_solution_time(std::numeric_limits<double>::quiet_NaN()),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_hier_math_ops(NULL),
      d_hier_math_ops_external(false),
      d_enable_logging(false)
{
    // intentionally blank
    return;
} // GeneralOperator()

GeneralOperator::~GeneralOperator()
{
    deallocateOperatorState();
    return;
} // ~GeneralOperator()

const std::string&
GeneralOperator::getName() const
{
    return d_object_name;
} // getName

bool
GeneralOperator::getIsInitialized() const
{
    return d_is_initialized;
} // getIsInitialized

void
GeneralOperator::setHomogeneousBc(bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

bool
GeneralOperator::getHomogeneousBc() const
{
    return d_homogeneous_bc;
} // getHomogeneousBc

void
GeneralOperator::setSolutionTime(double solution_time)
{
    d_solution_time = solution_time;
    return;
} // setSolutionTime

double
GeneralOperator::getSolutionTime() const
{
    return d_solution_time;
} // getSolutionTime

void
GeneralOperator::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
} // setTimeInterval

std::pair<double, double>
GeneralOperator::getTimeInterval() const
{
    return std::make_pair(d_current_time, d_new_time);
} // getTimeInterval

double
GeneralOperator::getDt() const
{
    return d_new_time - d_current_time;
} // getDt

void
GeneralOperator::setHierarchyMathOps(Pointer<HierarchyMathOps> hier_math_ops)
{
    d_hier_math_ops = hier_math_ops;
    d_hier_math_ops_external = d_hier_math_ops;
    return;
} // setHierarchyMathOps

Pointer<HierarchyMathOps>
GeneralOperator::getHierarchyMathOps() const
{
    return d_hier_math_ops;
} // getHierarchyMathOps

void
GeneralOperator::applyAdd(SAMRAIVectorReal<NDIM, double>& x,
                          SAMRAIVectorReal<NDIM, double>& y,
                          SAMRAIVectorReal<NDIM, double>& z)
{
    // Guard against the case that y == z.
    Pointer<SAMRAIVectorReal<NDIM, double> > zz = z.cloneVector(z.getName());
    zz->allocateVectorData();
    zz->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&z, false));
    apply(x, *zz);
    z.add(Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false), zz);
    zz->freeVectorComponents();
    return;
} // applyAdd

void
GeneralOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& /*in*/,
                                         const SAMRAIVectorReal<NDIM, double>& /*out*/)
{
    d_is_initialized = true;
    return;
} // initializeOperatorState

void
GeneralOperator::deallocateOperatorState()
{
    d_is_initialized = false;
    return;
} // deallocateOperatorState

void
GeneralOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& /*y*/)
{
    // intentionally blank
    return;
} // modifyRhsForBcs

void
GeneralOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& /*u*/)
{
    // intentionally blank
    return;
} // imposeSolBcs

void
GeneralOperator::setLoggingEnabled(bool enable_logging)
{
    d_enable_logging = enable_logging;
    return;
} // setLoggingEnabled

bool
GeneralOperator::getLoggingEnabled() const
{
    return d_enable_logging;
} // getLoggingEnabled

void
GeneralOperator::printClassData(std::ostream& stream)
{
    stream << "\n"
           << "object_name = " << d_object_name << "\n"
           << "is_initialized = " << d_is_initialized << "\n"
           << "homogeneous_bc = " << d_homogeneous_bc << "\n"
           << "solution_time = " << d_solution_time << "\n"
           << "current_time = " << d_current_time << "\n"
           << "new_time = " << d_new_time << "\n"
           << "hier_math_ops = " << d_hier_math_ops.getPointer() << "\n"
           << "hier_math_ops_external = " << d_hier_math_ops_external << "\n"
           << "enable_logging = " << d_enable_logging << "\n";
    return;
} // printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
