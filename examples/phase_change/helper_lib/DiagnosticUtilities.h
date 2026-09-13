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

#ifndef included_PhaseChangeExamples_DiagnosticUtilities
#define included_PhaseChangeExamples_DiagnosticUtilities

#include <ibamr/config.h>

#include <ibamr/PhaseChangeHierarchyIntegrator.h>

#include <HierarchyCellDataOpsReal.h>
#include <PatchHierarchy.h>

#include <fstream>
#include <initializer_list>
#include <string>

namespace PhaseChangeExamples
{
void open_diagnostic_file(std::ofstream& stream, const std::string& filename);
void allocate_diagnostic_data(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                              std::initializer_list<int> indices,
                              int coarsest_level,
                              int finest_level,
                              double time);
void deallocate_diagnostic_data(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                std::initializer_list<int> indices,
                                int coarsest_level,
                                int finest_level);

// Keep the diagnostic's original level range fixed through a run. Cases that
// regrid explicitly request scratch allocation again before using the product.
class PhaseMassDiagnostic
{
public:
    PhaseMassDiagnostic(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                        SAMRAI::tbox::Pointer<IBAMR::PhaseChangeHierarchyIntegrator> integrator,
                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> density,
                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> heaviside,
                        double time);
    ~PhaseMassDiagnostic();
    PhaseMassDiagnostic(const PhaseMassDiagnostic&) = delete;
    PhaseMassDiagnostic& operator=(const PhaseMassDiagnostic&) = delete;
    int getFinestLevel() const;
    int getHeavisideIndex() const;
    void allocateData(double time);
    void multiply();
    double integral(int weight_index) const;

private:
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double>> d_data_ops;
    int d_density_idx, d_heaviside_idx, d_mass_idx, d_finest_level;
};
} // namespace PhaseChangeExamples
#endif
