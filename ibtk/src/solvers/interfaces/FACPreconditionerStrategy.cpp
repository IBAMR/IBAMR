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

#include "ibtk/FACPreconditionerStrategy.h"

#include "Box.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "tbox/ConstPointer.h"

#include <ostream>
#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FACPreconditionerStrategy::FACPreconditionerStrategy(std::string object_name, bool homogeneous_bc)
    : d_object_name(std::move(object_name)), d_homogeneous_bc(homogeneous_bc)
{
    // intentionally blank
    return;
} // FACPreconditionerStrategy

const std::string&
FACPreconditionerStrategy::getName() const
{
    return d_object_name;
} // getName

bool
FACPreconditionerStrategy::getIsInitialized() const
{
    return d_is_initialized;
} // getIsInitialized

void
FACPreconditionerStrategy::setFACPreconditioner(ConstPointer<FACPreconditioner> preconditioner)
{
    d_preconditioner = preconditioner;
    return;
} // setFACPreconditioner

void
FACPreconditionerStrategy::setHomogeneousBc(bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

bool
FACPreconditionerStrategy::getHomogeneousBc() const
{
    return d_homogeneous_bc;
} // getHomogeneousBc

void
FACPreconditionerStrategy::setSolutionTime(double solution_time)
{
    d_solution_time = solution_time;
    return;
} // setSolutionTime

double
FACPreconditionerStrategy::getSolutionTime() const
{
    return d_solution_time;
} // getSolutionTime

void
FACPreconditionerStrategy::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
} // setTimeInterval

std::pair<double, double>
FACPreconditionerStrategy::getTimeInterval() const
{
    return std::make_pair(d_current_time, d_new_time);
} // getTimeInterval

double
FACPreconditionerStrategy::getDt() const
{
    return d_new_time - d_current_time;
} // getDt

void
FACPreconditionerStrategy::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& /*solution*/,
                                                   const SAMRAIVectorReal<NDIM, double>& /*rhs*/)
{
    d_is_initialized = true;
    return;
} // initializeOperatorState

void
FACPreconditionerStrategy::deallocateOperatorState()
{
    d_is_initialized = false;
    return;
} // deallocateOperatorState

void
FACPreconditionerStrategy::allocateScratchData()
{
    // intentionally blank
    return;
}

void
FACPreconditionerStrategy::deallocateScratchData()
{
    // intentionally blank
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

Pointer<SAMRAIVectorReal<NDIM, double> >
FACPreconditionerStrategy::getLevelSAMRAIVectorReal(const SAMRAIVectorReal<NDIM, double>& vec, int level_num) const
{
    Pointer<SAMRAIVectorReal<NDIM, double> > level_vec = new SAMRAIVectorReal<NDIM, double>(
        vec.getName() + "::level_" + std::to_string(level_num), vec.getPatchHierarchy(), level_num, level_num);
    for (int comp = 0; comp < vec.getNumberOfComponents(); ++comp)
    {
        level_vec->addComponent(
            vec.getComponentVariable(comp), vec.getComponentDescriptorIndex(comp), vec.getControlVolumeIndex(comp));
    }
    return level_vec;
} // getLevelSAMRAIVectorReal

void
FACPreconditionerStrategy::printClassData(std::ostream& stream)
{
    stream << "\n"
           << "object_name = " << d_object_name << "\n"
           << "is_initialized = " << d_is_initialized << "\n"
           << "homogeneous_bc = " << d_homogeneous_bc << "\n"
           << "solution_time = " << d_solution_time << "\n"
           << "current_time = " << d_current_time << "\n"
           << "new_time = " << d_new_time << "\n";
    return;
} // printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
