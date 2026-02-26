// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_samrai_compatibility_layer
#define included_IBTK_samrai_compatibility_layer

#include <ibtk/config.h>

#include <ibtk/samrai_compatibility_names.h>

#include <samrai_compatibility/tbox/SAMRAIConstPointer.h>
#include <samrai_compatibility/tbox/SAMRAIPointer.h>

#include <SAMRAIConstPointer.h>
#include <SAMRAIPointer.h>

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
using Pointer = SAMRAIPointer<T>;

template <class T>
using ConstPointer = SAMRAIConstPointer<T>;
} // namespace SAMRAICompatibility
} // namespace IBTK

#endif // #ifndef included_IBTK_samrai_compatibility_layer
