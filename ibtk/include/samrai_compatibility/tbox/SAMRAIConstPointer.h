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

#ifndef included_IBTK_samrai_compatibility_tbox_SAMRAIConstPointer
#define included_IBTK_samrai_compatibility_tbox_SAMRAIConstPointer

#include <samrai_compatibility/samrai_compatibility_detect.h>
#include <samrai_compatibility/tbox/SAMRAIPointer.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/tbox/ConstPointer.h>)
#include <SAMRAI/tbox/ConstPointer.h>
#define IBTK_SAMRAI_COMPAT_HAS_LEGACY_CONST_POINTER 1
#elif IBTK_SAMRAI_HAS_INCLUDE(<tbox/ConstPointer.h>)
#include <tbox/ConstPointer.h>
#define IBTK_SAMRAI_COMPAT_HAS_LEGACY_CONST_POINTER 1
#else
#define IBTK_SAMRAI_COMPAT_HAS_LEGACY_CONST_POINTER 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
#if (IBTK_SAMRAI_COMPAT_HAS_LEGACY_CONST_POINTER)
using SAMRAIConstPointer = SAMRAI::tbox::ConstPointer<T>;
#else
using SAMRAIConstPointer = SAMRAIPointer<const T>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_HAS_LEGACY_CONST_POINTER

#endif // #ifndef included_IBTK_samrai_compatibility_tbox_SAMRAIConstPointer
