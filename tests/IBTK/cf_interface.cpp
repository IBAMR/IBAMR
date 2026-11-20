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

#include <tbox/Logger.h>

#include <petscsys.h>

#include <ArrayData.h>
#include <BergerRigoutsos.h>
#include <BoundaryBox.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <RobinBcCoefStrategy.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ibtk/app_namespaces.h>

namespace
{
// Do not print file names or line numbers of an error to the output file.
class TestAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    void logMessage(const std::string& message, const std::string&, const int) override
    {
        // SAMRAI appends an extra NUL to the message string.
        SAMRAI::tbox::pout << "Program abort called with error message\n\n    " << message.c_str() << std::endl
                           << std::flush;
    }
};

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
                    x[d] = location_index % 2 == 0 ? x_lower[d] : x_upper[d];
                else
                    x[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
            }
            if (acoef_data) (*acoef_data)(i, 0) = 1.0;
            if (bcoef_data) (*bcoef_data)(i, 0) = 0.0;
            if (gcoef_data) (*gcoef_data)(i, 0) = d_fcn(x);
        }
    }

    IntVector<NDIM> numberOfExtensionsFillable() const override
    {
        return IntVector<NDIM>(128);
    }

private:
    std::function<double(const VectorNd&)> d_fcn;
};

// Value of the ghost cells before the fill, so that cells that nothing sets are equal on all patches.
constexpr double SENTINEL = -12345.678;

// The fine level: the union F of its boxes, and the domain, in the fine index space.
struct FineLevel
{
    BoxArray<NDIM> boxes;
    Box<NDIM> domain;
    Box<NDIM> coarse_domain;
    IntVector<NDIM> ratio;            // To level 0.
    IntVector<NDIM> ratio_to_coarser; // To the next coarser level.
    bool coarse_data_extended;        // Whether the fill sets the coarse data beyond the physical boundary.
};

// Cell data have axis -1. Side data have the faces normal to axis (axis >= 0), indexed by the cell above the face.
Box<NDIM>
index_box(const Box<NDIM>& cell_box, const int axis)
{
    Box<NDIM> box = cell_box;
    if (axis >= 0) box.upper(axis) += 1;
    return box;
}

bool
in_fine_level(const FineLevel& fine, const int axis, const hier::Index<NDIM>& i)
{
    for (int k = 0; k < fine.boxes.getNumberOfBoxes(); ++k)
    {
        if (index_box(fine.boxes[k], axis).contains(i)) return true;
    }
    return false;
}

int
floor_div(const int i, const int r)
{
    return i >= 0 ? i / r : -((-i + r - 1) / r);
}

// Cells (and faces) that the test checks: all of them if the fill extends the coarse data beyond the physical
// boundary, and otherwise those whose coarse stencil is far from the boundary, where nothing sets boundary values.
bool
away_from_boundary(const FineLevel& fine, const hier::Index<NDIM>& i)
{
    if (fine.coarse_data_extended) return true;
    for (int d = 0; d < NDIM; ++d)
    {
        const int c = floor_div(i(d), fine.ratio(d));
        if (c < fine.coarse_domain.lower(d) + 2 || c > fine.coarse_domain.upper(d) - 2) return false;
    }
    return true;
}

// The directions (d, s) along which a ghost cell (or face) i gets a coarse-fine value: i is not in F, the first cell
// (face) of F along the direction is at the distance k, at most the depth limit, and the cell (face) after it is also
// in F.
std::vector<std::pair<int, int>>
coarse_fine_terms(const FineLevel& fine, const int axis, const hier::Index<NDIM>& i, const int depth_limit)
{
    std::vector<std::pair<int, int>> terms;
    if (!index_box(fine.domain, axis).contains(i) || in_fine_level(fine, axis, i)) return terms;
    for (int d = 0; d < NDIM; ++d)
    {
        for (const int s : { -1, +1 })
        {
            const int limit = std::min(depth_limit, fine.ratio_to_coarser(d));
            for (int k = 1; k <= limit; ++k)
            {
                hier::Index<NDIM> j1 = i, j2 = i;
                j1(d) += k * s;
                j2(d) += (k + 1) * s;
                if (!in_fine_level(fine, axis, j1)) continue;
                if (in_fine_level(fine, axis, j2)) terms.emplace_back(d, s);
                break;
            }
        }
    }
    return terms;
}

// A ghost value of one patch.
struct GhostValue
{
    int axis;
    int idx[NDIM];
    int first_layer; // In the first ghost layer of this patch, within one cell of its tangential range.
    double value;
};

bool
less_key(const GhostValue& a, const GhostValue& b)
{
    if (a.axis != b.axis) return a.axis < b.axis;
    for (int d = NDIM - 1; d >= 0; --d)
    {
        if (a.idx[d] != b.idx[d]) return a.idx[d] < b.idx[d];
    }
    return false;
}

bool
same_key(const GhostValue& a, const GhostValue& b)
{
    return !less_key(a, b) && !less_key(b, a);
}

bool
same_bits(const double a, const double b)
{
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

std::string
value_string(const double v)
{
    std::ostringstream os;
    os.precision(17);
    os << v;
    return os.str();
}

std::string
key_string(const GhostValue& v)
{
    std::string s = "axis " + std::to_string(v.axis) + " index (";
    for (int d = 0; d < NDIM; ++d) s += (d ? "," : "") + std::to_string(v.idx[d]);
    return s + ")";
}

// Value of a patch datum at index i, for cell data (axis = -1) or side data.
double
get_value(Pointer<Patch<NDIM>> patch, const int data_idx, const int axis, const hier::Index<NDIM>& i)
{
    if (axis < 0)
    {
        Pointer<CellData<NDIM, double>> data = patch->getPatchData(data_idx);
        return (*data)(CellIndex<NDIM>(i));
    }
    Pointer<SideData<NDIM, double>> data = patch->getPatchData(data_idx);
    return (*data)(SideIndex<NDIM>(i, axis, SideIndex<NDIM>::Lower));
}

// Box of the cell indices (axis = -1) or the side indices in the ghost box of the patch data.
Box<NDIM>
ghost_index_box(Pointer<Patch<NDIM>> patch, const int data_idx, const int axis)
{
    Pointer<PatchData<NDIM>> data = patch->getPatchData(data_idx);
    return index_box(data->getGhostBox(), axis);
}

template <typename F>
void
for_each_index(const Box<NDIM>& box, const int axis, F f)
{
    for (BoxIterator<NDIM> bi(box); bi; bi++)
    {
        f(axis, bi());
    }
}
} // namespace

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
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
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        const int ghost_width = input_db->getIntegerWithDefault("GHOST_WIDTH", 1);
        // Without CF_FIXUP, the fill is the refine operator only, and only the interior values are checked.
        const bool cf_fixup = input_db->getBoolWithDefault("CF_FIXUP", true);
        // The fill sets the layers 1 to min(g, R) of the ghost cells.
        const int depth_limit = ghost_width;
        const std::string refine_op_name = input_db->getStringWithDefault("REFINE_OP", "NONE");
        const bool is_cell = input_db->getString("CENTERING").compare("CELL") == 0;
        // With EXTRAP_TYPE, the test fills the ghost cells with linear refinement and the extrapolation at the physical
        // boundary instead of the coarse-fine interpolation, for a linear function that both reproduce.
        const std::string extrap_type = input_db->getStringWithDefault("EXTRAP_TYPE", "NONE");
        const bool test_extrap = extrap_type != "NONE";
        // With BC_TYPE = "ROBIN", the component has Dirichlet data equal to the function. With UNRELATED_ROBIN_DATA, it
        // fills a linear function with Robin data unrelated to it, as the pressure does, which only the linear
        // extension of the coarse data reproduces.
        const bool use_robin = !test_extrap && input_db->getStringWithDefault("BC_TYPE", "EXTRAPOLATION") == "ROBIN";
        const bool unrelated_robin_data = input_db->getBoolWithDefault("UNRELATED_ROBIN_DATA", false);
        const bool extend_linearly = input_db->getBoolWithDefault("EXTEND_COARSE_DATA_LINEARLY", false);
        const bool Q_is_linear = test_extrap || unrelated_robin_data;

        // Create cell-centered or side-centered quantities, and initialize them with a quadratic function. The second
        // quantity is filled with a no-corners fill pattern.
        Pointer<Variable<NDIM>> Q_var, Q_nc_var;
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

        // Along a periodic direction, the function is a quadratic of the distance in the period, so that the faces on
        // the two sides of a periodic boundary have the same value.
        const IntVector<NDIM> periodic_shift = grid_geometry->getPeriodicShift(IntVector<NDIM>(1));
        auto Q_fcn = [&](const VectorNd& x) -> double
        {
            double val = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                if (periodic_shift(d) != 0)
                {
                    const double s = (x[d] - grid_geometry->getXLower()[d]) /
                                     (grid_geometry->getXUpper()[d] - grid_geometry->getXLower()[d]);
                    val += s * (1.0 - s);
                }
                else
                {
                    val += x[d] + (Q_is_linear ? 0.0 : x[d] * x[d]);
                    for (int e = d + 1; e < NDIM; ++e)
                    {
                        if (periodic_shift(e) == 0) val += Q_is_linear ? 0.0 : x[d] * x[e];
                    }
                }
            }
            return val;
        };

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int Q_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("CTX"), IntVector<NDIM>(ghost_width));
        const int Q_nc_idx =
            var_db->registerVariableAndContext(Q_nc_var, var_db->getContext("CTX"), IntVector<NDIM>(ghost_width));

        // With RESET_TO_OTHER_DESTINATION, the test builds the schedules to fill a copy of Q in place, then resets the
        // component to fill Q from the copy.
        const bool reset_to_other_destination = input_db->getBoolWithDefault("RESET_TO_OTHER_DESTINATION", false);
        const int Q_copy_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("COPY"), IntVector<NDIM>(ghost_width));
        // Initialize the AMR patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        if (input_db->keyExists("FINE_BOXES"))
        {
            // The boxes of the fine levels, spread over the processes. The refinement ratios are one integer for every
            // axis or one integer for each axis.
            auto read_ratio = [&](const std::string& key)
            {
                const tbox::Array<int> values = input_db->getIntegerArray(key);
                IntVector<NDIM> ratio;
                for (int d = 0; d < NDIM; ++d) ratio(d) = values[values.size() == 1 ? 0 : d];
                return ratio;
            };
            IntVector<NDIM> ratio_to_level_0(1);
            for (int ln = 1; ln <= 2; ++ln)
            {
                const std::string suffix = ln == 1 ? "" : "_2";
                if (!input_db->keyExists("FINE_BOXES" + suffix)) break;
                const tbox::Array<tbox::DatabaseBox> db_boxes = input_db->getDatabaseBoxArray("FINE_BOXES" + suffix);
                BoxArray<NDIM> fine_boxes(db_boxes.size());
                ProcessorMapping mapping(db_boxes.size());
                for (int k = 0; k < db_boxes.size(); ++k)
                {
                    fine_boxes[k] = Box<NDIM>(db_boxes[k]);
                    mapping.setProcessorAssignment(k, k % IBTK_MPI::getNodes());
                }
                ratio_to_level_0 *= read_ratio("REF_RATIO" + suffix);
                patch_hierarchy->makeNewPatchLevel(ln, ratio_to_level_0, fine_boxes, mapping);
            }
        }
        else
        {
            int tag_buffer = 1;
            int level_number = 0;
            bool done = false;
            while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
            {
                gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
                done = !patch_hierarchy->finerLevelExists(level_number);
                ++level_number;
            }
        }

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
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                for (const int data_idx : { Q_idx, Q_nc_idx })
                {
                    Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(data_idx);
                    Pointer<SideData<NDIM, double>> Q_sc_data = patch->getPatchData(data_idx);
                    if (Q_cc_data)
                    {
                        Q_cc_data->fillAll(SENTINEL);
                        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                        {
                            const CellIndex<NDIM>& idx = ci();
                            VectorNd x;
                            for (int d = 0; d < NDIM; ++d)
                                x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                            (*Q_cc_data)(idx) = Q_fcn(x);
                        }
                    }
                    else if (Q_sc_data)
                    {
                        Q_sc_data->fillAll(SENTINEL);
                        for (int axis = 0; axis < NDIM; ++axis)
                        {
                            for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                            {
                                const SideIndex<NDIM>& idx = si();
                                VectorNd x;
                                for (int d = 0; d < NDIM; ++d)
                                    x[d] = xlow[d] +
                                           dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
                                (*Q_sc_data)(idx) = Q_fcn(x);
                            }
                        }
                    }
                }
            }
        }

        // Fill the ghost cells with the coarse-fine interface operator, with the default fill pattern and with a
        // no-corners fill pattern.
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_comps(1);
        std::vector<std::unique_ptr<FunctionDirichletBc>> bcs;
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs;
        if (use_robin)
        {
            for (int c = 0; c < (is_cell ? 1 : NDIM); ++c)
            {
                bcs.push_back(std::make_unique<FunctionDirichletBc>(
                    unrelated_robin_data ?
                        std::function<double(const VectorNd&)>([](const VectorNd&) { return 1000.0; }) :
                        std::function<double(const VectorNd&)>(Q_fcn)));
                bc_coefs.push_back(bcs.back().get());
            }
        }
        ghost_cell_comps[0] = test_extrap ? InterpolationTransactionComponent(Q_idx,
                                                                              "CONSERVATIVE_LINEAR_REFINE",
                                                                              false,
                                                                              "CONSERVATIVE_COARSEN",
                                                                              extrap_type,
                                                                              false,
                                                                              nullptr /*bdry_conds*/) :
                                            InterpolationTransactionComponent(Q_idx,
                                                                              Q_idx,
                                                                              refine_op_name,
                                                                              cf_fixup,
                                                                              "NONE",
                                                                              "LINEAR",
                                                                              false,
                                                                              bc_coefs,
                                                                              nullptr /*fill_pattern*/);
        ghost_cell_comps[0].d_extend_coarse_data_linearly = extend_linearly;
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
                    if (Q_cc_data) Q_cc_data->fillAll(std::numeric_limits<double>::quiet_NaN());
                    if (Q_sc_data) Q_sc_data->fillAll(std::numeric_limits<double>::quiet_NaN());
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
            no_corners_pattern = new CellNoCornersFillPattern(ghost_width, /*overwrite_interior*/ false);
        else
            no_corners_pattern = new SideNoCornersFillPattern(ghost_width, /*overwrite_interior*/ false);
        std::vector<InterpolationTransactionComponent> no_corners_comps(1);
        no_corners_comps[0] = InterpolationTransactionComponent(
            Q_nc_idx, Q_nc_idx, refine_op_name, cf_fixup, "NONE", "LINEAR", false, bc_coefs, no_corners_pattern);
        no_corners_comps[0].d_extend_coarse_data_linearly = extend_linearly;
        HierarchyGhostCellInterpolation hier_ghost_cell_no_corners;
        hier_ghost_cell_no_corners.initializeOperatorState(no_corners_comps, patch_hierarchy);
        hier_ghost_cell_no_corners.fillData(0.0);

        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
        FineLevel fine;
        fine.boxes = level->getBoxes();
        fine.ratio = level->getRatio();
        fine.ratio_to_coarser = level->getRatioToCoarserLevel();
        fine.domain = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], level->getRatio());
        fine.coarse_domain = grid_geometry->getPhysicalDomain()[0];
        fine.coarse_data_extended = use_robin || ghost_width >= 2;
        std::vector<int> axes;
        if (is_cell)
        {
            axes.push_back(-1);
        }
        else
        {
            for (int d = 0; d < NDIM; ++d) axes.push_back(d);
        }

        // The exact value of the quadratic function at a ghost cell or face of a patch.
        auto exact_value = [&](Pointer<Patch<NDIM>> patch, const int axis, const hier::Index<NDIM>& i)
        {
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            const hier::Index<NDIM>& idx_low = patch->getBox().lower();
            VectorNd x;
            for (int d = 0; d < NDIM; ++d)
                x[d] = xlow[d] + dx[d] * (static_cast<double>(i(d) - idx_low(d)) + (d == axis ? 0.0 : 0.5));
            return Q_fcn(x);
        };

        std::ostringstream failures;

        // Check that every ghost cell (or face) that is not in the fine level and that has fine-level cells (or faces)
        // on one side at distances one and two gets the exact value of the quadratic function. Such a cell is at the
        // coarse-fine interface of some patch of the level, whatever its position relative to the patch whose ghost box
        // contains it: beyond the end of the patch, deep in its ghost region, or at a convex or a concave corner.
        if (cf_fixup && !test_extrap)
        {
            int num_checked = 0, num_wrong = 0;
            std::string wrong_cells;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int axis : axes)
                {
                    for_each_index(
                        ghost_index_box(patch, Q_idx, axis),
                        axis,
                        [&](const int ax, const hier::Index<NDIM>& i)
                        {
                            if (!away_from_boundary(fine, i) || coarse_fine_terms(fine, ax, i, depth_limit).empty())
                                return;
                            ++num_checked;
                            const double value = get_value(patch, Q_idx, ax, i);
                            if (!IBTK::abs_equal_eps(value, exact_value(patch, ax, i)))
                            {
                                if (num_wrong++ < 60)
                                {
                                    std::ostringstream os;
                                    os << "    axis " << ax << " index " << i << " patch " << patch->getBox()
                                       << " computed " << value << " exact " << exact_value(patch, ax, i) << "\n";
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
        if (cf_fixup && !test_extrap)
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
                                       if (interior_box.contains(i) || !index_box(fine.domain, ax).contains(i) ||
                                           !away_from_boundary(fine, i) || in_fine_level(fine, ax, i))
                                           return;
                                       GhostValue v;
                                       v.axis = ax;
                                       for (int d = 0; d < NDIM; ++d) v.idx[d] = i(d);
                                       v.first_layer = ring_box.contains(i) ? 1 : 0;
                                       v.value = get_value(patch, Q_idx, ax, i);
                                       mine.push_back(v);
                                   });
                }
            }
            // Gather the values of all patches on all processes.
            const int nproc = IBTK_MPI::getNodes();
            const int nbytes = static_cast<int>(mine.size() * sizeof(GhostValue));
            std::vector<int> counts(nproc), displs(nproc);
            MPI_Allgather(&nbytes, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
            int total = 0;
            for (int r = 0; r < nproc; ++r)
            {
                displs[r] = total;
                total += counts[r];
            }
            std::vector<GhostValue> all(total / sizeof(GhostValue));
            MPI_Allgatherv(
                mine.data(), nbytes, MPI_BYTE, all.data(), counts.data(), displs.data(), MPI_BYTE, MPI_COMM_WORLD);
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
                        if (all[last].first_layer) in_set = true;
                        ++last;
                    }
                    if (last - first > 1 && (set == 'A' || in_set))
                    {
                        ++num_shared;
                        bool differ = false;
                        for (std::size_t k = first + 1; k < last; ++k)
                        {
                            if (!same_bits(all[first].value, all[k].value)) differ = true;
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
        }

        // Check that every patch that contains a cell or a face of the fine level, as an interior value or as a ghost
        // value, has the same value there, to the bit. A face on the boundary of a patch along a coarse-fine interface
        // belongs to the fine level, so the ghost copies of it in other patches, on this or another process, must equal
        // the value of the patch that contains it, also for patches that do not touch that patch. Cells and faces at
        // periodic images are the same cell or face.
        if (!test_extrap)
        {
            const IntVector<NDIM> shift = grid_geometry->getPeriodicShift(level->getRatio());
            std::vector<GhostValue> mine;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int axis : axes)
                {
                    for_each_index(ghost_index_box(patch, Q_idx, axis),
                                   axis,
                                   [&](const int ax, const hier::Index<NDIM>& i)
                                   {
                                       hier::Index<NDIM> j = i;
                                       for (int d = 0; d < NDIM; ++d)
                                       {
                                           if (shift(d) == 0) continue;
                                           const int lo = fine.domain.lower(d);
                                           j(d) = lo + ((i(d) - lo) % shift(d) + shift(d)) % shift(d);
                                       }
                                       if (!in_fine_level(fine, ax, j)) return;
                                       GhostValue v;
                                       v.axis = ax;
                                       for (int d = 0; d < NDIM; ++d) v.idx[d] = j(d);
                                       v.first_layer = index_box(patch->getBox(), ax).contains(i) ? 0 : 1;
                                       v.value = get_value(patch, Q_idx, ax, i);
                                       mine.push_back(v);
                                   });
                }
            }
            const int nproc = IBTK_MPI::getNodes();
            const int nbytes = static_cast<int>(mine.size() * sizeof(GhostValue));
            std::vector<int> counts(nproc), displs(nproc);
            MPI_Allgather(&nbytes, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
            int total = 0;
            for (int r = 0; r < nproc; ++r)
            {
                displs[r] = total;
                total += counts[r];
            }
            std::vector<GhostValue> all(total / sizeof(GhostValue));
            MPI_Allgatherv(
                mine.data(), nbytes, MPI_BYTE, all.data(), counts.data(), displs.data(), MPI_BYTE, MPI_COMM_WORLD);
            std::stable_sort(all.begin(), all.end(), less_key);

            int num_shared = 0, num_differ = 0;
            std::string examples;
            for (std::size_t first = 0; first < all.size();)
            {
                std::size_t last = first;
                while (last < all.size() && same_key(all[first], all[last])) ++last;
                if (last - first > 1)
                {
                    ++num_shared;
                    bool differ = false;
                    for (std::size_t k = first + 1; k < last; ++k)
                    {
                        if (!same_bits(all[first].value, all[k].value)) differ = true;
                    }
                    if (differ)
                    {
                        if (num_differ < 40)
                        {
                            examples += "    " + key_string(all[first]) + ":";
                            for (std::size_t k = first; k < last; ++k)
                            {
                                examples += " " + value_string(all[k].value) + (all[k].first_layer ? "*" : "");
                            }
                            examples += "\n";
                        }
                        ++num_differ;
                    }
                }
                first = last;
            }
            pout << "Compared " << num_shared << (is_cell ? " cells" : " faces")
                 << " of the fine level that several patches contain: " << num_differ << " differ.\n";
            if (num_differ > 0)
            {
                failures << "  patches disagree on " << num_differ << (is_cell ? " cells" : " faces")
                         << " of the fine level (* marks a ghost copy):\n"
                         << examples;
            }
        }

        // Check that the cells that the no-corners fill fills have the values of the default fill, to the bit: the
        // cells of the face slabs of each patch, as wide as the ghost width.
        if (cf_fixup && !test_extrap)
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
                                               if (!index_box(fine.domain, ax).contains(i) ||
                                                   !away_from_boundary(fine, i))
                                                   return;
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

        // Check that the fill does not change the values in the interiors of the patches, including the faces of their
        // boundaries, to the bit.
        if (!test_extrap)
        {
            int num_compared = 0, num_changed = 0;
            std::string first_changed;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int data_idx : { Q_idx, Q_nc_idx })
                {
                    for (const int axis : axes)
                    {
                        for_each_index(index_box(patch->getBox(), axis),
                                       axis,
                                       [&](const int ax, const hier::Index<NDIM>& i)
                                       {
                                           ++num_compared;
                                           const double value = get_value(patch, data_idx, ax, i);
                                           if (!same_bits(value, exact_value(patch, ax, i)) && num_changed++ == 0)
                                           {
                                               std::ostringstream os;
                                               os << "axis " << ax << " index " << i << " patch data " << data_idx
                                                  << " value " << value_string(value);
                                               first_changed = os.str();
                                           }
                                       });
                    }
                }
            }
            IBTK_MPI::sumReduction(&num_compared, 1);
            int num_changed_all = num_changed;
            IBTK_MPI::sumReduction(&num_changed_all, 1);
            pout << "Compared " << num_compared
                 << " interior values of the patches with their values before the fill: " << num_changed_all
                 << " changed.\n";
            if (num_changed_all > 0)
            {
                failures << "  " << num_changed_all
                         << " interior values changed by the fill; first on process with a change: " << first_changed
                         << "\n";
            }
        }

        // Check the ghost cells between the lower boundary and a patch within the ghost cell width of it. The cells are
        // filled from a temporary coarse patch whose ghost cells are set by extrapolation.
        if (test_extrap)
        {
            int num_checked = 0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                const double* const xlow = pgeom->getXLower();
                const hier::Index<NDIM>& idx_low = patch->getBox().lower();
                Pointer<CellData<NDIM, double>> Q_cc_data = patch->getPatchData(Q_idx);
                TBOX_ASSERT(Q_cc_data);
                // The ghost cells below the patch, away from its corners and from the other boundaries.
                Box<NDIM> below = patch->getBox();
                below.upper(1) = below.lower(1) - 1;
                below.lower(1) -= ghost_width;
                for (int d = 0; d < NDIM; ++d)
                {
                    if (d != 1)
                    {
                        below.lower(d) += 2;
                        below.upper(d) -= 2;
                    }
                }
                below = below * fine.domain;
                if (below.empty()) continue;
                for (CellIterator<NDIM> ci(below); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                    ++num_checked;
                    if (!IBTK::abs_equal_eps((*Q_cc_data)(idx)-Q_fcn(x), 0.0))
                    {
                        failures << "  incorrect ghost value on cell " << idx << ", patch " << patch->getBox()
                                 << ": computed " << (*Q_cc_data)(idx) << ", exact " << Q_fcn(x) << "\n";
                    }
                }
            }
            pout << "Checked " << num_checked << " values below patches.\n";
        }

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
