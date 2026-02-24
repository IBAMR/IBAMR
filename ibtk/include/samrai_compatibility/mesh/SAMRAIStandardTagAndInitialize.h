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

#ifndef included_IBTK_samrai_compatibility_mesh_SAMRAIStandardTagAndInitialize
#define included_IBTK_samrai_compatibility_mesh_SAMRAIStandardTagAndInitialize

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

// SAMRAI may provide helper functions in headers/sources that are intentionally
// unused in some build configurations.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#endif

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/mesh/StandardTagAndInitialize.h")
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_StandardTagAndInitialize 1
#else
#include <StandardTagAndInitialize.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_StandardTagAndInitialize 0
#endif

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_StandardTagAndInitialize)
using SAMRAIStandardTagAndInitialize = SAMRAI::mesh::StandardTagAndInitialize;
#else
using SAMRAIStandardTagAndInitialize = SAMRAI::mesh::StandardTagAndInitialize<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_StandardTagAndInitialize

#endif // #ifndef included_IBTK_samrai_compatibility_mesh_SAMRAIStandardTagAndInitialize
