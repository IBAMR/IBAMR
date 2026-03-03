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

#ifndef included_IBTK_samrai_compatibility_detect
#define included_IBTK_samrai_compatibility_detect

#if defined(__has_include)
#define IBTK_SAMRAI_HAS_INCLUDE(header) __has_include(header)
#else
#define IBTK_SAMRAI_HAS_INCLUDE(header) 0
#endif

#endif // #ifndef included_IBTK_samrai_compatibility_detect
