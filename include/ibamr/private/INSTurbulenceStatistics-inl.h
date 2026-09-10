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

#ifndef included_IBAMR_INSTurbulenceStatistics_inl
#define included_IBAMR_INSTurbulenceStatistics_inl

#include <ibamr/config.h>

#include <ibamr/INSTurbulenceStatistics.h>

namespace IBAMR
{
inline INSTurbulenceStatistics::INSTurbulenceStatistics(std::string object_name, const double statistics_start_time)
    : d_object_name(std::move(object_name)), d_statistics_start_time(statistics_start_time)
{
    // intentionally blank
}

inline bool
INSTurbulenceStatistics::shouldUpdateStatistics(const double data_time)
{
    return data_time >= d_statistics_start_time;
}

inline void
INSTurbulenceStatistics::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM>> /*visit_writer*/)
{
    // intentionally blank
}

inline void
INSTurbulenceStatistics::setupPlotData(const double /*data_time*/,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> /*hierarchy*/,
                                       SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> /*hier_math_ops*/)
{
    // intentionally blank
}

inline void
INSTurbulenceStatistics::deallocatePlotData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> /*hierarchy*/)
{
    // intentionally blank
}
} // namespace IBAMR

#endif
