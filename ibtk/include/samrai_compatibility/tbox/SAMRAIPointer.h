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

#ifndef included_IBTK_samrai_compatibility_tbox_SAMRAIPointer
#define included_IBTK_samrai_compatibility_tbox_SAMRAIPointer

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/tbox/Pointer.h")
#include <SAMRAI/tbox/Pointer.h>
#define IBTK_SAMRAI_COMPAT_HAS_LEGACY_POINTER 1
#elif IBTK_SAMRAI_HAS_INCLUDE("tbox/Pointer.h")
#include <tbox/Pointer.h>
#define IBTK_SAMRAI_COMPAT_HAS_LEGACY_POINTER 1
#else
#include <memory>
#define IBTK_SAMRAI_COMPAT_HAS_LEGACY_POINTER 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
#if (IBTK_SAMRAI_COMPAT_HAS_LEGACY_POINTER)
using SAMRAIPointer = SAMRAI::tbox::Pointer<T>;
#else
using SAMRAIPointer = std::shared_ptr<T>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_HAS_LEGACY_POINTER

#endif // #ifndef included_IBTK_samrai_compatibility_tbox_SAMRAIPointer
