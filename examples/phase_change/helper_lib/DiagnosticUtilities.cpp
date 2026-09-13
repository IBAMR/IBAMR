// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <VariableDatabase.h>

#include "DiagnosticUtilities.h"

#include <ibamr/app_namespaces.h>

namespace PhaseChangeExamples
{
void
open_diagnostic_file(std::ofstream& stream, const std::string& filename)
{
    stream.open(filename, std::ios_base::out | std::ios_base::app);
    stream.precision(16);
    stream.setf(std::ios::fixed, std::ios::floatfield);
}

void
allocate_diagnostic_data(Pointer<PatchHierarchy<NDIM>> hierarchy,
                         const std::initializer_list<int> indices,
                         const int coarsest_level,
                         const int finest_level,
                         const double time)
{
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const int idx : indices)
        {
            if (!level->checkAllocated(idx))
            {
                level->allocatePatchData(idx, time);
            }
        }
    }
}

void
deallocate_diagnostic_data(Pointer<PatchHierarchy<NDIM>> hierarchy,
                           const std::initializer_list<int> indices,
                           const int coarsest_level,
                           const int finest_level)
{
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        for (const int idx : indices)
        {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(idx);
        }
    }
    for (const int idx : indices)
    {
        VariableDatabase<NDIM>::getDatabase()->removePatchDataIndex(idx);
    }
}

PhaseMassDiagnostic::PhaseMassDiagnostic(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                         Pointer<PhaseChangeHierarchyIntegrator> integrator,
                                         Pointer<CellVariable<NDIM, double>> density,
                                         Pointer<CellVariable<NDIM, double>> heaviside,
                                         const double time)
    : d_hierarchy(hierarchy)
{
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    d_density_idx = var_db->mapVariableAndContextToIndex(density, integrator->getCurrentContext());
    d_heaviside_idx = var_db->mapVariableAndContextToIndex(heaviside, integrator->getCurrentContext());
    d_mass_idx = var_db->registerClonedPatchDataIndex(heaviside, d_heaviside_idx);
    d_finest_level = hierarchy->getFinestLevelNumber();
    allocateData(time);
    d_data_ops = new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, 0, d_finest_level);
}

PhaseMassDiagnostic::~PhaseMassDiagnostic()
{
    deallocate_diagnostic_data(d_hierarchy, { d_mass_idx }, 0, d_finest_level);
}

int
PhaseMassDiagnostic::getFinestLevel() const
{
    return d_finest_level;
}

int
PhaseMassDiagnostic::getHeavisideIndex() const
{
    return d_heaviside_idx;
}

void
PhaseMassDiagnostic::allocateData(const double time)
{
    allocate_diagnostic_data(d_hierarchy, { d_mass_idx }, 0, d_finest_level, time);
}

void
PhaseMassDiagnostic::multiply()
{
    d_data_ops->multiply(d_mass_idx, d_density_idx, d_heaviside_idx);
}

double
PhaseMassDiagnostic::integral(const int weight_index) const
{
    return d_data_ops->integral(d_mass_idx, weight_index);
}
} // namespace PhaseChangeExamples
