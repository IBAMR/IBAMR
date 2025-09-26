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

template <class T>
FACPreconditionerStrategy<T>::FACPreconditionerStrategy(std::string object_name, bool homogeneous_bc)
    : d_object_name(std::move(object_name)), d_homogeneous_bc(homogeneous_bc)
{
    // intentionally blank
    return;
} // FACPreconditionerStrategy

template <class T>
const std::string&
FACPreconditionerStrategy<T>::getName() const
{
    return d_object_name;
} // getName

template <class T>
bool
FACPreconditionerStrategy<T>::getIsInitialized() const
{
    return d_is_initialized;
} // getIsInitialized

template <class T>
void
FACPreconditionerStrategy<T>::setFACPreconditioner(ConstPointer<FACPreconditioner<T> > preconditioner)
{
    d_preconditioner = preconditioner;
    return;
} // setFACPreconditioner

template <class T>
void
FACPreconditionerStrategy<T>::setHomogeneousBc(bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

template <class T>
bool
FACPreconditionerStrategy<T>::getHomogeneousBc() const
{
    return d_homogeneous_bc;
} // getHomogeneousBc

template <class T>
void
FACPreconditionerStrategy<T>::setSolutionTime(double solution_time)
{
    d_solution_time = solution_time;
    return;
} // setSolutionTime

template <class T>
double
FACPreconditionerStrategy<T>::getSolutionTime() const
{
    return d_solution_time;
} // getSolutionTime

template <class T>
void
FACPreconditionerStrategy<T>::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
} // setTimeInterval

template <class T>
std::pair<double, double>
FACPreconditionerStrategy<T>::getTimeInterval() const
{
    return std::make_pair(d_current_time, d_new_time);
} // getTimeInterval

template <class T>
double
FACPreconditionerStrategy<T>::getDt() const
{
    return d_new_time - d_current_time;
} // getDt

template <class T>
void
FACPreconditionerStrategy<T>::initializeOperatorState(const SAMRAIVectorReal<NDIM, T>& /*solution*/,
                                                      const SAMRAIVectorReal<NDIM, T>& /*rhs*/)
{
    d_is_initialized = true;
    return;
} // initializeOperatorState

template <class T>
void
FACPreconditionerStrategy<T>::deallocateOperatorState()
{
    d_is_initialized = false;
    return;
} // deallocateOperatorState

template <class T>
void
FACPreconditionerStrategy<T>::allocateScratchData()
{
    // intentionally blank
    return;
}

template <class T>
void
FACPreconditionerStrategy<T>::deallocateScratchData()
{
    // intentionally blank
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

template <class T>
Pointer<SAMRAIVectorReal<NDIM, T> >
FACPreconditionerStrategy<T>::getLevelSAMRAIVectorReal(const SAMRAIVectorReal<NDIM, T>& vec, int level_num) const
{
    Pointer<SAMRAIVectorReal<NDIM, T> > level_vec = new SAMRAIVectorReal<NDIM, T>(
        vec.getName() + "::level_" + std::to_string(level_num), vec.getPatchHierarchy(), level_num, level_num);
    for (int comp = 0; comp < vec.getNumberOfComponents(); ++comp)
    {
        level_vec->addComponent(
            vec.getComponentVariable(comp), vec.getComponentDescriptorIndex(comp), vec.getControlVolumeIndex(comp));
    }
    return level_vec;
} // getLevelSAMRAIVectorReal

template <class T>
void
FACPreconditionerStrategy<T>::printClassData(std::ostream& stream)
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

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class FACPreconditionerStrategy<float>;
template class FACPreconditionerStrategy<double>;

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
