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

#ifndef included_PhaseChangeTestUtilities
#define included_PhaseChangeTestUtilities

#include <ibamr/config.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/IBTK_MPI.h>

#include <CellVariable.h>
#include <SideData.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

namespace IBAMR
{
class PhaseChangeHierarchyIntegrator;
}

class RefinementRegion
{
public:
    /*! \brief Move the refinement band every timestep of length dt. */
    explicit RefinementRegion(double dt);
    /*! \brief Return the center of the band at the supplied time. */
    double getCenter(double time) const;
    /*! \brief Exchange the two band positions. */
    void swapRegions();

private:
    double d_dt;
    bool d_swap_regions = false;
};

// Count actual invalidations of the phase-change operators when the mesh changes.
template <class Integrator>
class RegridCountingIntegrator : public Integrator
{
public:
    using Integrator::Integrator;
    /*! \brief Return the number of hierarchy reconfigurations. */
    int getConfigurationResetCount() const;

    /*! \brief Return the number of changes to the local fine-level boxes. */
    int getMeshChangeCount() const;

protected:
    /*! \brief Record reconfigurations after resetting the inherited operators. */
    void resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                                int coarsest_level,
                                                int finest_level) override;

private:
    int d_configuration_reset_count = 0;
    int d_mesh_change_count = 0;
    std::string d_fine_boxes;
};

void tag_moving_refinement_region(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM>> hierarchy,
                                  int level_number,
                                  double time,
                                  int tag_idx,
                                  bool initial_time,
                                  bool uses_richardson_extrapolation,
                                  void* ctx);

// Compare every cell and side value, including covered coarse cells, against
// the uninterrupted trajectory written by attest's first run of this input.
void check_restart_fields(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                          const std::vector<int>& cell_indices,
                          const std::vector<int>& side_indices,
                          int step,
                          double time,
                          bool from_restart,
                          std::ostream& results);

void check_liquid_fraction_tags(SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy);

void check_liquid_fraction_gradient(SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> fraction_var,
                                    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> gradient_var,
                                    std::ostream& results);

void check_divergence_source_transfer(SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> root_integrator,
                                      SAMRAI::tbox::Pointer<IBAMR::PhaseChangeHierarchyIntegrator> phase_integrator,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                      RefinementRegion& region,
                                      std::ostream& results);

#include "PhaseChangeTestUtilities-inl.h"

#endif
