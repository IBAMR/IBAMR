// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/private/CachedGhostCellFill.h>
#include <ibtk/private/InterfaceFaceGhostCopy.h>

#include <tbox/Logger.h>

#include <petscsys.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <RobinBcCoefStrategy.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <new>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../coarse_fine_ghost_values.h"
#include "../tests.h"

#include <ibtk/app_namespaces.h>

using namespace coarse_fine_ghost_values;

// Fill the ghost cells (or faces) of a quadratic function at coarse-fine interfaces with the quadratic coarse-fine
// interpolation, with the default fill pattern and with a no-corners fill pattern, and check on the finest level that
//
//   - every ghost value that is at a coarse-fine interface of some patch equals the quadratic function;
//   - all patches, on all processes, agree to the bit on the value of every ghost value outside the fine level;
//   - the no-corners fill gives the values of the default fill in the face slabs of each patch;
//   - the interior values are unchanged; and
//   - the ghost values equal the values of the Reference block of the input to 1e-12. The values are sums and maxima of
//     the ghost values, over each ghost location once, so that they do not depend on the partition of the fine level
//     into patches, on the processes that own the patches, or on the order of the summation.
//
// With REPLACE_LEVEL_0 = TRUE, the test instead fills levels 1 and 2 with a cached fill operator, replaces level 0 by a
// new PatchLevel object with new data, and checks that the cached fill then gives the ghost values of a fresh fill.
//
// With REPLACE_BC_OBJECT = TRUE, the test instead fills with a cached fill operator and one Robin coefficient object,
// destroys the object, constructs a different object with different boundary data in the same storage, and checks that
// the next cached fill, which uses its operator, gives the ghost values of a fresh fill with the new object.
//
// With RESET_FACE_COPY_PATTERN = TRUE and side-centered data, the test instead releases the fill pattern of an
// InterfaceFaceGhostCopy, which has a pattern of the kind OTHER, and resets the object with a different pattern that is
// constructed in the storage of the released one. The reset must stop with an error. The test returns normally if it
// does not.
//
// With RESET_ADDS = "EXTRAPOLATION" or "ROBIN" and cell-centered data, the test instead sets up a fill with neither a
// physical boundary extrapolation nor Robin coefficients, and resets its component with the given one. The reset must
// stop with an error. The test returns normally if it does not.
//
// For each run, the input filename must be given on the command line. In all cases, the command line is:
//
//    executable <input file name>

namespace
{
// Dirichlet data equal to a function at the physical boundary, evaluated at the cell centers along the boundary. For
// side-centered data, the boundary operator shifts the patch geometry.
class FunctionDirichletBc : public RobinBcCoefStrategy<NDIM>
{
public:
    explicit FunctionDirichletBc(std::function<double(const VectorNd&)> fcn) : d_fcn(std::move(fcn))
    {
    }

    void setBcCoefs(Pointer<ArrayData<NDIM, double>>& acoef_data,
                    Pointer<ArrayData<NDIM, double>>& bcoef_data,
                    Pointer<ArrayData<NDIM, double>>& gcoef_data,
                    const Pointer<SAMRAI::hier::Variable<NDIM>>& /*variable*/,
                    const Patch<NDIM>& patch,
                    const BoundaryBox<NDIM>& bdry_box,
                    double /*fill_time*/) const override
    {
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const x_upper = pgeom->getXUpper();
        const double* const dx = pgeom->getDx();
        const hier::Index<NDIM>& patch_lower = patch.getBox().lower();
        const unsigned int location_index = bdry_box.getLocationIndex();
        const unsigned int normal_axis = location_index / 2;
        const Box<NDIM>& box = (acoef_data ? acoef_data->getBox() : gcoef_data->getBox());
        for (Box<NDIM>::Iterator b(box); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            VectorNd x;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d == normal_axis)
                {
                    x[d] = location_index % 2 == 0 ? x_lower[d] : x_upper[d];
                }
                else
                {
                    x[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                }
            }
            if (acoef_data)
            {
                (*acoef_data)(i, 0) = 1.0;
            }
            if (bcoef_data)
            {
                (*bcoef_data)(i, 0) = 0.0;
            }
            if (gcoef_data)
            {
                (*gcoef_data)(i, 0) = d_fcn(x);
            }
        }
    }

    IntVector<NDIM> numberOfExtensionsFillable() const override
    {
        return IntVector<NDIM>(128);
    }

private:
    std::function<double(const VectorNd&)> d_fcn;
};

// Summary of the ghost values at the distinct ghost locations outside the fine level that are away from the physical
// boundary. The ghost locations are cells (or faces) in the ghost region of some patch. The coarse-fine locations are
// those that get a coarse-fine value.
struct GhostSummary
{
    double num_locations = 0.0;
    double num_coarse_fine = 0.0;
    double sum = 0.0;
    double sum_squares = 0.0;
    double max_abs = 0.0;
    double weighted_sum = 0.0;         // Weight 1 + sum_d (d + 1) (i(d) - lower(d)) detects a swap of two values.
    double weighted_sum_squares = 0.0; // Sets the size of the terms of the weighted sum.
};

// Print a quantity with enough digits to distinguish partitions, and add a failure to failures if it differs from the
// value of the same key in the Reference block by more than 1e-12, relative to the largest of the reference value, 1,
// and term_size. The roundoff of a sum whose terms cancel scales with the root of the sum of the squares of its terms.
// The caller passes that root as term_size.
void
check_against_reference(Pointer<Database> reference_db,
                        const std::string& key,
                        const std::string& label,
                        const double value,
                        const double term_size,
                        std::ostringstream& failures)
{
    std::ostringstream os;
    os << std::setprecision(14) << label << " = " << value << '\n';
    if (!reference_db->keyExists(key))
    {
        os << "  the Reference block has no key " << key << '\n';
        failures << "  the Reference block of the input has no key " << key << "\n";
    }
    else
    {
        const double reference = reference_db->getDouble(key);
        if (!(std::abs(value - reference) <= 1.0e-12 * std::max({ 1.0, std::abs(reference), term_size })))
        {
            os << "  differs from the reference value " << reference << '\n';
            failures << "  " << label << " = " << std::setprecision(14) << value << " differs from the reference value "
                     << reference << "\n";
        }
    }
    pout << os.str();
} // check_against_reference

void
check_summary_against_reference(Pointer<Database> input_db, const GhostSummary& s, std::ostringstream& failures)
{
    if (!input_db->keyExists("Reference"))
    {
        failures << "  the input has no Reference block\n";
        return;
    }
    Pointer<Database> db = input_db->getDatabase("Reference");
    check_against_reference(db, "ghost_locations", "ghost locations", s.num_locations, 0.0, failures);
    check_against_reference(
        db, "coarse_fine_locations", "coarse-fine ghost locations", s.num_coarse_fine, 0.0, failures);
    check_against_reference(
        db, "ghost_sum", "sum of coarse-fine ghost values", s.sum, std::sqrt(s.sum_squares), failures);
    check_against_reference(
        db, "ghost_sum_squares", "sum of squares of coarse-fine ghost values", s.sum_squares, 0.0, failures);
    check_against_reference(db, "ghost_max", "largest absolute coarse-fine ghost value", s.max_abs, 0.0, failures);
    check_against_reference(db,
                            "ghost_weighted_sum",
                            "weighted sum of coarse-fine ghost values",
                            s.weighted_sum,
                            std::sqrt(s.weighted_sum_squares),
                            failures);
} // check_summary_against_reference

// Fill levels 1 and 2 with a cached fill operator, replace level 0 by a new PatchLevel object with the same boxes and
// different data, and fill again. The cached fill must read the new level 0, as a freshly built fill does. Only cell
// data are supported.
void
check_cached_fill_after_replacing_level_0(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                          const int Q_idx,
                                          const int Q_fresh_idx,
                                          const std::string& refine_op_name,
                                          const std::function<double(const VectorNd&)>& Q_fcn)
{
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    const auto make_comps = [&](const int idx)
    {
        return std::vector<InterpolationTransactionComponent>(
            1, InterpolationTransactionComponent(idx, refine_op_name, true, "NONE", "LINEAR"));
    };
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    if (finest_ln != 2)
    {
        TBOX_ERROR("cf_interface: REPLACE_LEVEL_0 needs three levels\n");
    }

    IBTK::detail::CachedGhostCellFill cache;
    cache.get(make_comps(Q_idx), patch_hierarchy, 1, finest_ln, false).fillData(0.0);

    // Replace level 0, with the same boxes, mapping, and ratio, and give all levels new data.
    Pointer<PatchLevel<NDIM>> old_level = patch_hierarchy->getPatchLevel(0);
    patch_hierarchy->makeNewPatchLevel(
        0, old_level->getRatio(), old_level->getBoxes(), old_level->getProcessorMapping());
    patch_hierarchy->getPatchLevel(0)->allocatePatchData(Q_idx);
    patch_hierarchy->getPatchLevel(0)->allocatePatchData(Q_fresh_idx);
    fill_quadratic(patch_hierarchy, { Q_idx }, [&](const VectorNd& x) { return Q_fcn(x) + 1.0; });
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            patch->getPatchData(Q_fresh_idx)->copy(*patch->getPatchData(Q_idx));
        }
    }

    cache.get(make_comps(Q_idx), patch_hierarchy, 1, finest_ln, false).fillData(0.0);
    HierarchyGhostCellInterpolation fresh;
    fresh.initializeOperatorState(make_comps(Q_fresh_idx), patch_hierarchy, 1, finest_ln);
    fresh.fillData(0.0);

    int num_filled = 0;
    double max_diff = 0.0;
    for (int ln = 1; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> cached_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double>> fresh_data = patch->getPatchData(Q_fresh_idx);
            for (CellIterator<NDIM> ci(cached_data->getGhostBox()); ci; ci++)
            {
                if (patch->getBox().contains(ci()))
                {
                    continue;
                }
                const double value = (*fresh_data)(ci(), 0);
                num_filled += value != SENTINEL;
                max_diff = std::max(max_diff, std::abs((*cached_data)(ci(), 0) - value));
            }
        }
    }
    IBTK_MPI::sumReduction(&num_filled, 1);
    max_diff = IBTK_MPI::maxReduction(max_diff);
    pout << "Compared " << num_filled
         << " ghost cells of levels 1 and 2 that a fresh fill sets: the largest difference "
         << "between the cached fill and the fresh fill is " << max_diff << ".\n";
    if (max_diff != 0.0)
    {
        TBOX_ERROR("cf_interface: the cached fill differs from a fresh fill after level 0 was replaced\n");
    }
} // check_cached_fill_after_replacing_level_0

// Fill with a cached fill operator and a Robin coefficient object, destroy the object, construct a different object at
// the same address, and fill again. The cached fill must use the new object, as a freshly built fill does, and reuse
// its operator. Q_old_idx receives the ghost values of the first fill and Q_fresh_idx those of the fresh fill. Only
// cell data are supported.
void
check_cached_fill_after_replacing_bc_object(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                            const int Q_idx,
                                            const int Q_old_idx,
                                            const int Q_fresh_idx,
                                            const std::string& refine_op_name,
                                            const std::function<double(const VectorNd&)>& Q_fcn)
{
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    const auto make_comps = [&](const int idx, RobinBcCoefStrategy<NDIM>* const bc_coef)
    {
        return std::vector<InterpolationTransactionComponent>(
            1, InterpolationTransactionComponent(idx, refine_op_name, true, "NONE", "LINEAR", false, bc_coef));
    };
    const auto copy_data = [&](const int dst_idx, const int src_idx)
    {
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                patch->getPatchData(dst_idx)->copy(*patch->getPatchData(src_idx));
            }
        }
    };
    copy_data(Q_fresh_idx, Q_idx);

    alignas(FunctionDirichletBc) unsigned char storage[sizeof(FunctionDirichletBc)];
    auto* const bc = new (storage) FunctionDirichletBc(Q_fcn);
    IBTK::detail::CachedGhostCellFill cache;
    HierarchyGhostCellInterpolation& fill = cache.get(make_comps(Q_idx, bc), patch_hierarchy);
    fill.fillData(0.0);
    copy_data(Q_old_idx, Q_idx);

    // The new object has the same address and different data, so the Dirichlet values and the ghost values at the
    // physical boundary differ by 2 from those of the first fill.
    bc->~FunctionDirichletBc();
    auto* new_bc = new (storage) FunctionDirichletBc([&](const VectorNd& x) { return Q_fcn(x) + 1.0; });
    const bool same_operator = &cache.get(make_comps(Q_idx, new_bc), patch_hierarchy) == &fill;
    fill.fillData(0.0);

    HierarchyGhostCellInterpolation fresh;
    fresh.initializeOperatorState(make_comps(Q_fresh_idx, new_bc), patch_hierarchy);
    fresh.fillData(0.0);

    int num_changed = 0;
    double max_diff = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> cached_data = patch->getPatchData(Q_idx);
            Pointer<CellData<NDIM, double>> old_data = patch->getPatchData(Q_old_idx);
            Pointer<CellData<NDIM, double>> fresh_data = patch->getPatchData(Q_fresh_idx);
            for (CellIterator<NDIM> ci(cached_data->getGhostBox()); ci; ci++)
            {
                if (patch->getBox().contains(ci()))
                {
                    continue;
                }
                num_changed += (*cached_data)(ci(), 0) != (*old_data)(ci(), 0);
                max_diff = std::max(max_diff, std::abs((*cached_data)(ci(), 0) - (*fresh_data)(ci(), 0)));
            }
        }
    }
    IBTK_MPI::sumReduction(&num_changed, 1);
    max_diff = IBTK_MPI::maxReduction(max_diff);
    new_bc->~FunctionDirichletBc();
    pout << "The second fill changed " << num_changed << " ghost cells: the largest difference between the cached fill "
         << "and the fresh fill is " << max_diff << ".\n";
    if (!same_operator)
    {
        TBOX_ERROR("cf_interface: the cached fill built a new operator after the Robin coefficient object changed\n");
    }
    if (num_changed == 0)
    {
        TBOX_ERROR("cf_interface: the new Robin coefficient object changed no ghost value\n");
    }
    if (max_diff != 0.0)
    {
        TBOX_ERROR(
            "cf_interface: the cached fill differs from a fresh fill after the Robin coefficient object was "
            "replaced\n");
    }
} // check_cached_fill_after_replacing_bc_object

// A fill pattern of the kind OTHER whose objects use two fixed slots, so that a new object takes the address of a
// released one if the slot is free.
class SlotFillPattern : public SideNoCornersFillPattern
{
public:
    using SideNoCornersFillPattern::SideNoCornersFillPattern;

    static void* operator new(std::size_t size)
    {
        for (int k = 0; k < 2; ++k)
        {
            if (!s_used[k] && size <= sizeof(s_slots[k]))
            {
                s_used[k] = true;
                return s_slots[k];
            }
        }
        throw std::bad_alloc();
    }

    static void operator delete(void* ptr)
    {
        for (int k = 0; k < 2; ++k)
        {
            if (ptr == s_slots[k])
            {
                s_used[k] = false;
            }
        }
    }

private:
    inline static bool s_used[2] = { false, false };
    alignas(std::max_align_t) inline static unsigned char s_slots[2][256];
};

// Construct an InterfaceFaceGhostCopy with a fill pattern, release the pattern, and reset the object with a different
// pattern that is constructed in the storage of the released one. The reset must stop, because the pattern differs,
// although the addresses are equal. If the reset returns, this function returns normally.
void
check_face_copy_reset_with_new_pattern(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                       const int Q_idx,
                                       const int ghost_width)
{
    Pointer<VariableFillPattern<NDIM>> pattern = new SlotFillPattern(ghost_width, false);
    IBTK::detail::InterfaceFaceGhostCopy face_copy(
        patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber(), Q_idx, Q_idx, pattern, false);
    pattern.setNull();
    Pointer<VariableFillPattern<NDIM>> new_pattern = new SlotFillPattern(ghost_width, false);
    face_copy.resetPatchDataIndices(Q_idx, Q_idx, new_pattern);
} // check_face_copy_reset_with_new_pattern

// Set up a fill with neither a physical boundary extrapolation nor Robin coefficients, then reset its component with a
// LINEAR extrapolation (kind EXTRAPOLATION) or with the Robin coefficient object bc_coef (kind ROBIN). The reset must
// stop. If the reset returns, this function returns normally.
void
check_reset_adding_physical_boundary_fill(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                          const int Q_idx,
                                          const std::string& refine_op_name,
                                          const std::string& kind,
                                          RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> comps(
        1, InterpolationTransactionComponent(Q_idx, refine_op_name, true, "NONE", "NONE"));
    HierarchyGhostCellInterpolation hier_ghost_cell;
    hier_ghost_cell.initializeOperatorState(comps, patch_hierarchy);
    if (kind == "EXTRAPOLATION")
    {
        comps[0].d_phys_bdry_extrap_type = "LINEAR";
    }
    else
    {
        comps[0].d_robin_bc_coefs = { bc_coef };
    }
    hier_ghost_cell.resetTransactionComponents(comps);
} // check_reset_adding_physical_boundary_fill
} // namespace

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Suppress the SAMRAI warning for patches smaller than the ghost cell width, because the warning contains a path.
    Logger::getInstance()->setWarning(false);

    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ghost_cells.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        const int ghost_width = input_db->getIntegerWithDefault("GHOST_WIDTH", 1);
        // The fill sets the layers 1 to min(g, R) of the ghost cells.
        const int depth_limit = ghost_width;
        const std::string refine_op_name = input_db->getStringWithDefault("REFINE_OP", "NONE");
        const bool is_cell = input_db->getString("CENTERING").compare("CELL") == 0;
        // With BC_TYPE = "ROBIN", the component has Dirichlet data equal to the function. With UNRELATED_ROBIN_DATA, it
        // fills a linear function with Robin data unrelated to it, as the pressure does, which only the linear
        // extension of the coarse data reproduces.
        const bool use_robin = input_db->getStringWithDefault("BC_TYPE", "EXTRAPOLATION") == "ROBIN";
        const bool unrelated_robin_data = input_db->getBoolWithDefault("UNRELATED_ROBIN_DATA", false);
        const bool extend_linearly = input_db->getBoolWithDefault("EXTEND_COARSE_DATA_LINEARLY", false);

        // Create cell-centered or side-centered quantities, and initialize them with a quadratic function. The second
        // quantity is filled with a no-corners fill pattern.
        Pointer<hier::Variable<NDIM>> Q_var, Q_nc_var;
        if (is_cell)
        {
            Q_var = new CellVariable<NDIM, double>("Q");
            Q_nc_var = new CellVariable<NDIM, double>("Q_no_corners");
        }
        else
        {
            Q_var = new SideVariable<NDIM, double>("Q");
            Q_nc_var = new SideVariable<NDIM, double>("Q_no_corners");
        }

        const std::function<double(const VectorNd&)> linear_fcn = [](const VectorNd& x)
        {
            double val = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                val += x[d];
            }
            return val;
        };
        const std::function<double(const VectorNd&)> Q_fcn =
            unrelated_robin_data ? linear_fcn : make_quadratic(grid_geometry);

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int Q_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("CTX"), IntVector<NDIM>(ghost_width));
        const int Q_nc_idx =
            var_db->registerVariableAndContext(Q_nc_var, var_db->getContext("CTX"), IntVector<NDIM>(ghost_width));

        // With RESET_TO_OTHER_DESTINATION, the test builds the schedules to fill a copy of Q in place, then resets the
        // component to fill Q from the copy.
        const bool reset_to_other_destination = input_db->getBoolWithDefault("RESET_TO_OTHER_DESTINATION", false);
        // With CLONE_COPY, the copy is a cloned index, which has a variable and no context.
        const bool clone_copy = input_db->getBoolWithDefault("CLONE_COPY", false);
        const int Q_copy_idx =
            clone_copy ?
                var_db->registerClonedPatchDataIndex(Q_var, Q_idx) :
                var_db->registerVariableAndContext(Q_var, var_db->getContext("COPY"), IntVector<NDIM>(ghost_width));
        // Initialize the AMR patch hierarchy.
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = make_hierarchy(app_initializer, grid_geometry);

        // Allocate patch data
        int coarsest_ln = 0, finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_idx);
            level->allocatePatchData(Q_nc_idx);
            level->allocatePatchData(Q_copy_idx);
        }

        // Fill in initial data
        fill_quadratic(patch_hierarchy, { Q_idx, Q_nc_idx }, Q_fcn);

        if (input_db->getBoolWithDefault("REPLACE_LEVEL_0", false))
        {
            // Q_copy receives the result of the fresh fill.
            check_cached_fill_after_replacing_level_0(patch_hierarchy, Q_idx, Q_copy_idx, refine_op_name, Q_fcn);
            return 0;
        }
        if (input_db->getBoolWithDefault("REPLACE_BC_OBJECT", false))
        {
            // Q_copy receives the result of the first fill and Q_nc the result of the fresh fill.
            check_cached_fill_after_replacing_bc_object(
                patch_hierarchy, Q_idx, Q_copy_idx, Q_nc_idx, refine_op_name, Q_fcn);
            return 0;
        }

        if (input_db->getBoolWithDefault("RESET_FACE_COPY_PATTERN", false))
        {
            check_face_copy_reset_with_new_pattern(patch_hierarchy, Q_idx, ghost_width);
            return 0;
        }
        if (input_db->keyExists("RESET_ADDS"))
        {
            FunctionDirichletBc bc(Q_fcn);
            check_reset_adding_physical_boundary_fill(
                patch_hierarchy, Q_idx, refine_op_name, input_db->getString("RESET_ADDS"), &bc);
            return 0;
        }

        // Fill the ghost cells with the coarse-fine interface operator, with the default fill pattern and with a
        // no-corners fill pattern.
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<std::unique_ptr<FunctionDirichletBc>> bcs;
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs;
        if (use_robin)
        {
            const std::function<double(const VectorNd&)> bc_fcn =
                unrelated_robin_data ? std::function<double(const VectorNd&)>([](const VectorNd&) { return 1000.0; }) :
                                       Q_fcn;
            for (int c = 0; c < (is_cell ? 1 : NDIM); ++c)
            {
                bcs.push_back(std::make_unique<FunctionDirichletBc>(bc_fcn));
                bc_coefs.push_back(bcs.back().get());
            }
        }
        std::vector<InterpolationTransactionComponent> ghost_cell_comps(1);
        ghost_cell_comps[0] = InterpolationTransactionComponent(Q_idx,
                                                                Q_idx,
                                                                refine_op_name,
                                                                true,
                                                                "NONE",
                                                                "LINEAR",
                                                                false,
                                                                bc_coefs,
                                                                nullptr /*fill_pattern*/,
                                                                "LINEAR",
                                                                extend_linearly);
        HierarchyGhostCellInterpolation hier_ghost_cell;
        if (reset_to_other_destination)
        {
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    patch->getPatchData(Q_copy_idx)->copy(*patch->getPatchData(Q_idx));
                    Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(Q_idx);
                    Pointer<SideData<NDIM, double>> Q_sc_data = patch->getPatchData(Q_idx);
                    if (Q_cc_data)
                    {
                        Q_cc_data->fillAll(std::numeric_limits<double>::quiet_NaN());
                    }
                    if (Q_sc_data)
                    {
                        Q_sc_data->fillAll(std::numeric_limits<double>::quiet_NaN());
                    }
                }
            }
            ghost_cell_comps[0].d_dst_data_idx = Q_copy_idx;
            ghost_cell_comps[0].d_src_data_idx = Q_copy_idx;
            hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
            ghost_cell_comps[0].d_dst_data_idx = Q_idx;
            hier_ghost_cell.resetTransactionComponents(ghost_cell_comps);
        }
        else
        {
            hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        }
        hier_ghost_cell.fillData(0.0);

        Pointer<VariableFillPattern<NDIM>> no_corners_pattern;
        if (is_cell)
        {
            no_corners_pattern = new CellNoCornersFillPattern(ghost_width, /*overwrite_interior*/ false);
        }
        else
        {
            no_corners_pattern = new SideNoCornersFillPattern(ghost_width, /*overwrite_interior*/ false);
        }
        std::vector<InterpolationTransactionComponent> no_corners_comps(1);
        no_corners_comps[0] = InterpolationTransactionComponent(Q_nc_idx,
                                                                Q_nc_idx,
                                                                refine_op_name,
                                                                true,
                                                                "NONE",
                                                                "LINEAR",
                                                                false,
                                                                bc_coefs,
                                                                no_corners_pattern,
                                                                "LINEAR",
                                                                extend_linearly);
        HierarchyGhostCellInterpolation hier_ghost_cell_no_corners;
        hier_ghost_cell_no_corners.initializeOperatorState(no_corners_comps, patch_hierarchy);
        hier_ghost_cell_no_corners.fillData(0.0);

        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
        FineLevel fine = make_fine_level(patch_hierarchy, grid_geometry);
        fine.coarse_data_extended = use_robin || ghost_width >= 2;
        const std::vector<int> axes = data_axes(is_cell);

        // The exact value of the quadratic function at a ghost cell or face of a patch.
        auto exact = [&](Pointer<Patch<NDIM>> patch, const int axis, const hier::Index<NDIM>& i)
        { return exact_value(patch, axis, i, Q_fcn); };

        std::ostringstream failures;
        GhostSummary summary;

        // Check that every ghost cell (or face) that is not in the fine level and that has fine-level cells (or faces)
        // on one side at distances one and two gets the exact value of the quadratic function. Such a cell is at the
        // coarse-fine interface of some patch of the level, whatever its position relative to the patch whose ghost box
        // contains it: beyond the end of the patch, deep in its ghost region, or at a convex or a concave corner.
        {
            int num_checked = 0, num_wrong = 0;
            std::string wrong_cells;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int axis : axes)
                {
                    for_each_index(ghost_index_box(patch, Q_idx, axis),
                                   axis,
                                   [&](const int ax, const hier::Index<NDIM>& i)
                                   {
                                       if (!away_from_boundary(fine, i) ||
                                           coarse_fine_terms(fine, ax, i, depth_limit).empty())
                                       {
                                           return;
                                       }
                                       ++num_checked;
                                       const double value = get_value(patch, Q_idx, ax, i);
                                       if (!IBTK::abs_equal_eps(value, exact(patch, ax, i)))
                                       {
                                           if (num_wrong++ < 60)
                                           {
                                               std::ostringstream os;
                                               os << "    axis " << ax << " index " << i << " patch " << patch->getBox()
                                                  << " computed " << value << " exact " << exact(patch, ax, i) << "\n";
                                               wrong_cells += os.str();
                                           }
                                       }
                                   });
                }
            }
            IBTK_MPI::sumReduction(&num_checked, 1);
            int num_wrong_all = num_wrong;
            IBTK_MPI::sumReduction(&num_wrong_all, 1);
            pout << "Checked " << num_checked
                 << " ghost values at coarse-fine interfaces against the quadratic function.\n";
            if (num_wrong_all > 0)
            {
                failures << "  " << num_wrong_all
                         << " coarse-fine ghost values differ from the quadratic function (the first ones on this "
                            "process):\n"
                         << wrong_cells;
            }
        }

        // Check that all patches agree, to the bit, on the value of every ghost cell (or face) that is not in the fine
        // level. Set A is every such ghost cell. Set B is the cells that are in the first ghost layer of some patch,
        // within one cell of the tangential range of that patch; it is a diagnostic, a subset of set A.
        {
            std::vector<GhostValue> mine;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int axis : axes)
                {
                    const Box<NDIM> ring_box = index_box(Box<NDIM>::grow(patch->getBox(), IntVector<NDIM>(1)), axis);
                    const Box<NDIM> interior_box = index_box(patch->getBox(), axis);
                    for_each_index(ghost_index_box(patch, Q_idx, axis),
                                   axis,
                                   [&](const int ax, const hier::Index<NDIM>& i)
                                   {
                                       if (interior_box.contains(i) || !in_boxes(fine.domain, ax, i) ||
                                           !away_from_boundary(fine, i) || in_fine_level(fine, ax, i))
                                       {
                                           return;
                                       }
                                       GhostValue v;
                                       v.axis = ax;
                                       for (int d = 0; d < NDIM; ++d)
                                       {
                                           v.idx[d] = i(d);
                                       }
                                       v.first_layer = ring_box.contains(i) ? 1 : 0;
                                       v.value = get_value(patch, Q_idx, ax, i);
                                       mine.push_back(v);
                                   });
                }
            }
            // Gather the values of all patches on all processes.
            std::vector<GhostValue> all = gather_all(mine);
            std::stable_sort(all.begin(), all.end(), less_key);

            int num_differ_a = 0;
            std::string examples;
            for (const char set : { 'B', 'A' })
            {
                int num_shared = 0, num_differ = 0;
                std::string set_examples;
                for (std::size_t first = 0; first < all.size();)
                {
                    std::size_t last = first;
                    bool in_set = false;
                    while (last < all.size() && same_key(all[first], all[last]))
                    {
                        if (all[last].first_layer)
                        {
                            in_set = true;
                        }
                        ++last;
                    }
                    if (last - first > 1 && (set == 'A' || in_set))
                    {
                        ++num_shared;
                        bool differ = false;
                        for (std::size_t k = first + 1; k < last; ++k)
                        {
                            if (!same_bits(all[first].value, all[k].value))
                            {
                                differ = true;
                            }
                        }
                        if (differ)
                        {
                            if (num_differ < 40)
                            {
                                set_examples += "    " + key_string(all[first]) + ":";
                                for (std::size_t k = first; k < last; ++k)
                                {
                                    set_examples += " " + value_string(all[k].value) + (all[k].first_layer ? "*" : "");
                                }
                                set_examples += "\n";
                            }
                            ++num_differ;
                        }
                    }
                    first = last;
                }
                pout << "Compared " << num_shared << " ghost cells that several patches contain, set " << set
                     << (set == 'A' ? " (all ghost cells outside the fine level)" :
                                      " (first ghost layer, within one cell of the tangential range)")
                     << ": " << num_differ << " differ.\n";
                if (set == 'A')
                {
                    num_differ_a = num_differ;
                    examples = set_examples;
                }
            }
            if (num_differ_a > 0)
            {
                failures << "  patches disagree on " << num_differ_a
                         << " ghost values in set A (* marks the first ghost layer of the patch):\n"
                         << examples;
            }

            // Summarize the value at each ghost location once, in the order of the sorted keys, so that every
            // partition sums the same values in the same order. Every process has all the values.
            for (std::size_t k = 0; k < all.size(); ++k)
            {
                if (k > 0 && same_key(all[k - 1], all[k]))
                {
                    continue;
                }
                summary.num_locations += 1.0;
                hier::Index<NDIM> i;
                double weight = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    i(d) = all[k].idx[d];
                    weight += (d + 1) * (i(d) - fine.domain.getBoundingBox().lower(d));
                }
                if (coarse_fine_terms(fine, all[k].axis, i, depth_limit).empty())
                {
                    continue;
                }
                const double v = all[k].value;
                summary.num_coarse_fine += 1.0;
                summary.sum += v;
                summary.sum_squares += v * v;
                summary.max_abs = std::max(summary.max_abs, std::abs(v));
                summary.weighted_sum += weight * v;
                summary.weighted_sum_squares += weight * weight * v * v;
            }
        }

        check_fine_level_agreement(level, fine, grid_geometry, Q_idx, axes, is_cell, failures);

        // Check that the cells that the no-corners fill fills have the values of the default fill, to the bit: the
        // cells of the face slabs of each patch, as wide as the ghost width.
        {
            int num_compared = 0, num_differ = 0;
            std::string first_differ;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int axis : axes)
                {
                    const Box<NDIM> patch_index_box = index_box(patch->getBox(), axis);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        for (const int side : { -1, +1 })
                        {
                            Box<NDIM> slab = patch_index_box;
                            if (side < 0)
                            {
                                slab.lower(d) = patch_index_box.lower(d) - ghost_width;
                                slab.upper(d) = patch_index_box.lower(d) - 1;
                            }
                            else
                            {
                                slab.lower(d) = patch_index_box.upper(d) + 1;
                                slab.upper(d) = patch_index_box.upper(d) + ghost_width;
                            }
                            for_each_index(slab,
                                           axis,
                                           [&](const int ax, const hier::Index<NDIM>& i)
                                           {
                                               if (!in_boxes(fine.domain, ax, i) || !away_from_boundary(fine, i))
                                               {
                                                   return;
                                               }
                                               ++num_compared;
                                               const double v = get_value(patch, Q_idx, ax, i);
                                               const double v_nc = get_value(patch, Q_nc_idx, ax, i);
                                               if (!same_bits(v, v_nc) && num_differ++ == 0)
                                               {
                                                   std::ostringstream os;
                                                   os << "axis " << ax << " index " << i << " default fill " << v
                                                      << " no-corners fill " << v_nc;
                                                   first_differ = os.str();
                                               }
                                           });
                        }
                    }
                }
            }
            IBTK_MPI::sumReduction(&num_compared, 1);
            int num_differ_all = num_differ;
            IBTK_MPI::sumReduction(&num_differ_all, 1);
            pout << "Compared " << num_compared
                 << " ghost values of the no-corners fill with the default fill: " << num_differ_all << " differ.\n";
            if (num_differ_all > 0)
            {
                failures << "  " << num_differ_all
                         << " ghost values of the no-corners fill differ from the default fill; first on process "
                            "with a difference: "
                         << first_differ << "\n";
            }
        }

        check_interior_unchanged(level, { Q_idx, Q_nc_idx }, axes, Q_fcn, failures);

        check_summary_against_reference(input_db, summary, failures);

        if (!failures.str().empty())
        {
            TBOX_ERROR("cf_interface:\n" << failures.str());
        }

        // Deallocate patch data
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_idx);
            level->deallocatePatchData(Q_nc_idx);
            level->deallocatePatchData(Q_copy_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
