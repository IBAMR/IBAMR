// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/GeneralOperator.h"
#include "ibtk/HierarchyMathOps.h"

#include "Box.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include <ostream>
#include <string>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

GeneralOperator::GeneralOperator(std::string object_name, bool homogeneous_bc)
    : d_object_name(std::move(object_name)), d_homogeneous_bc(homogeneous_bc)
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
