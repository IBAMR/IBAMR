// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"

#include "CellVariable.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "tbox/Timer.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

// Local includes
#include "CartLaplaceOperator.h"
#include "KDTree.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;

static unsigned int num_pts = 10;
} // namespace

unsigned int CartLaplaceOperator::s_num_ghost_cells = 2;

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartLaplaceOperator::CartLaplaceOperator(const std::string& object_name,
                                         BoundaryMesh* bdry_mesh,
                                         const double dist_to_bdry)
    : PETScLinearAugmentedOperator(object_name, false), d_mesh(bdry_mesh), d_dist_to_bdry(dist_to_bdry)
{
    // Setup Timers.
    IBTK_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBTK::LaplaceOperator::apply()");
                 t_initialize_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::LaplaceOperator::initializeOperatorState()");
                 t_deallocate_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::LaplaceOperator::deallocateOperatorState()"););
    return;
} // LaplaceOperator()

CartLaplaceOperator::~CartLaplaceOperator()
{
    if (d_is_initialized)
    {
        this->deallocateOperatorState();
    }
    return;
} // ~LaplaceOperator()

void
CartLaplaceOperator::setupBeforeApply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<CellVariable<NDIM, double> > x_cc_var = x.getComponentVariable(comp);
        Pointer<CellVariable<NDIM, double> > y_cc_var = y.getComponentVariable(comp);
        if (!x_cc_var || !y_cc_var)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  encountered non-cell centered vector components" << std::endl);
        }
        Pointer<CellDataFactory<NDIM, double> > x_factory = x_cc_var->getPatchDataFactory();
        Pointer<CellDataFactory<NDIM, double> > y_factory = y_cc_var->getPatchDataFactory();
        TBOX_ASSERT(x_factory);
        TBOX_ASSERT(y_factory);
        const unsigned int x_depth = x_factory->getDefaultDepth();
        const unsigned int y_depth = y_factory->getDefaultDepth();
        TBOX_ASSERT(x_depth == y_depth);
        if (x_depth != d_bc_coefs.size() || y_depth != d_bc_coefs.size())
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  each vector component must have data depth == " << d_bc_coefs.size() << "\n"
                                     << "  since d_bc_coefs.size() == " << d_bc_coefs.size() << std::endl);
        }
    }
#endif

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent x_component(x.getComponentDescriptorIndex(comp),
                                                      DATA_REFINE_TYPE,
                                                      USE_CF_INTERPOLATION,
                                                      DATA_COARSEN_TYPE,
                                                      BDRY_EXTRAP_TYPE,
                                                      CONSISTENT_TYPE_2_BDRY,
                                                      d_bc_coefs,
                                                      d_fill_pattern);
        transaction_comps.push_back(x_component);
    }
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);
}

void
CartLaplaceOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBTK_TIMER_START(t_apply);

    setupBeforeApply(x, y);
    // Compute the action of the operator.
    for (unsigned int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();

            Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(d_ls_idx);
            for (int comp = 0; comp < d_ncomp; ++comp)
            {
                Pointer<CellData<NDIM, double> > x_data = patch->getPatchData(x.getComponentDescriptorIndex(comp));
                Pointer<CellData<NDIM, double> > y_data = patch->getPatchData(y.getComponentDescriptorIndex(comp));

                for (CellIterator<NDIM> ci(box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    const double ls_val = (*ls_data)(idx);
                    // Check if ls_val is within some specified tolerance of the boundary
                    if (ls_val < 0.0 && std::abs(ls_val) >= d_dist_to_bdry)
                    {
                        // We want to use finite differences.
                        double diff = 0.0;
                        for (int dir = 0; dir < NDIM; ++dir)
                        {
                            IntVector<NDIM> dirs(0);
                            dirs(dir) = 1;
                            diff += ((*x_data)(idx + dirs) - 2.0 * (*x_data)(idx) + (*x_data)(idx - dirs)) /
                                    (dx[dir] * dx[dir]);
                        }
                        (*y_data)(idx) = diff;
                    }
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_apply);
    return;
} // apply

void
CartLaplaceOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                             const SAMRAIVectorReal<NDIM, double>& out)
{
    IBTK_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in.cloneVector(in.getName());
    d_b = out.cloneVector(out.getName());

    // Setup operator state.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();

    d_ncomp = in.getNumberOfComponents();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
    TBOX_ASSERT(d_ncomp == out.getNumberOfComponents());
#endif

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops =
            new HierarchyMathOps(d_object_name + "::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_hier_math_ops);
#endif
    }

    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.clear();
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent component(d_x->getComponentDescriptorIndex(comp),
                                                    DATA_REFINE_TYPE,
                                                    USE_CF_INTERPOLATION,
                                                    DATA_COARSEN_TYPE,
                                                    BDRY_EXTRAP_TYPE,
                                                    CONSISTENT_TYPE_2_BDRY,
                                                    d_bc_coefs);
        d_transaction_comps.push_back(component);
    }

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy, d_coarsest_ln, d_finest_ln);

    // Sort Lag DOFs
    sortLagDOFsToCells();

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
CartLaplaceOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_operator_state);

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();
    d_fill_pattern.setNull();

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Delete the solution and rhs vectors.
    d_x->freeVectorComponents();
    d_x.setNull();

    d_b->freeVectorComponents();
    d_b.setNull();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

void
CartLaplaceOperator::printPtMap(std::ostream& os)
{
    const int ln = d_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
    unsigned int patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
    {
        os << "On patch number: " << patch_num << "\n";
        os << "There are " << d_base_pt_vec[patch_num].size() << " key-value pairs present\n";
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        for (size_t i = 0; i < d_base_pt_vec[patch_num].size(); ++i)
        {
            const UPoint& pt = d_base_pt_vec[patch_num][i];
            const std::vector<UPoint>& pt_vec = d_pair_pt_vec[patch_num][i];
            os << "  Looking at point:\n" << pt << "\n";
            os << "  Has points: \n";
            for (const auto& pt_from_vec : pt_vec) os << pt_from_vec << "\n";
        }
        os << "\n";
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////
void
CartLaplaceOperator::applyToLagDOFs(const int x_idx, const int y_idx)
{
    // Now apply the RBF reconstruction operator near the boundary mesh
    // Assume structure is on finest level
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(d_hierarchy->getFinestLevelNumber());
    unsigned int patch_num = 0;
    auto rbf = [](const double r) -> double { return r * r * r; };
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        // First loop through Cartesian grid cells.
        const Box<NDIM>& box = patch->getBox();
        Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(d_ls_idx);
        Pointer<CellData<NDIM, double> > x_data = patch->getPatchData(x_idx);
        Pointer<CellData<NDIM, double> > y_data = patch->getPatchData(y_idx);

        for (CellIterator<NDIM> ci(box); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            const double ls = (*ls_data)(idx);
            if (ls > 0.0 || std::abs(ls) > d_dist_to_bdry) continue;
            // We are on an index that needs RBF treatment.
            // We need to collect points near this one
            // Using LINEAR polynomials with r^3, so we need 5 points
            int poly_size = NDIM + 1;
        }
    }
}

void
CartLaplaceOperator::sortLagDOFsToCells()
{
    d_idx_node_vec.clear();
    d_idx_node_ghost_vec.clear();
    d_base_pt_vec.clear();
    d_pair_pt_vec.clear();
    // Assume structure is on finest level
    int ln = d_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
    int num_patches = level->getNumberOfPatches();
    d_idx_node_vec.resize(num_patches);
    d_idx_node_ghost_vec.resize(num_patches);
    d_base_pt_vec.resize(num_patches);
    d_pair_pt_vec.resize(num_patches);

    // Loop through nodes
    const MeshBase::node_iterator it_end = d_mesh->nodes_end();
    for (auto it = d_mesh->nodes_begin(); it != it_end; ++it)
    {
        Node* const node = *it;
        // Get CellIndex of the node
        VectorNd node_pt;
        for (unsigned int d = 0; d < NDIM; ++d) node_pt[d] = (*node)(d);
        const CellIndex<NDIM>& idx =
            IndexUtilities::getCellIndex(node_pt, d_hierarchy->getGridGeometry(), level->getRatio());
        // Sort nodes into patches
        unsigned int patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            if (box.contains(idx)) d_idx_node_vec[patch_num].push_back(node);
            Box<NDIM> ghost_box = box;
            ghost_box.grow(IntVector<NDIM>(s_num_ghost_cells));
            if (ghost_box.contains(idx)) d_idx_node_ghost_vec[patch_num].push_back(node);
        }
    }
    // At this point, each node is associated with a patch
    // We now find the nearest neighbor of every point that needs an RBF reconstruction
    // We do this by finding a KD tree
    unsigned int patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
    {
        if (d_idx_node_vec[patch_num].size() == 0) continue;
        // We are on a patch that has points. We need to form a KD tree.
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(d_ls_idx);

        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        // First start by collecting all points into a vector
        std::vector<UPoint> pts;
        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            if ((*ls_data)(idx) < 0.0) pts.push_back(UPoint(patch, idx));
        }
        // Now the points in d_idx_node_vec
        for (const auto& node : d_idx_node_ghost_vec[patch_num]) pts.push_back(UPoint(patch, node));

        // Now create KD tree
        tree::KDTree<UPoint> tree(pts);
        // We have a tree, now we need to find closest points for each point.
        // Start with Eulerian points
        size_t i = 0;
        for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            const double ls = (*ls_data)(idx);
            if (ls < 0.0 && std::abs(ls) < d_dist_to_bdry)
            {
                std::vector<int> idx_vec;
                std::vector<double> distance_vec;
                d_base_pt_vec[patch_num].push_back(UPoint(patch, idx));
                d_pair_pt_vec[patch_num].push_back({});
                tree.knnSearch(UPoint(patch, idx), num_pts, idx_vec, distance_vec);
                for (const auto& idx_in_pts : idx_vec) d_pair_pt_vec[patch_num][i].push_back(pts[idx_in_pts]);
                ++i;
            }
        }
        // Now do Lagrangian points
        for (const auto& node : d_idx_node_vec[patch_num])
        {
            std::vector<int> idx_vec;
            std::vector<double> distance_vec;
            d_base_pt_vec[patch_num].push_back(UPoint(patch, node));
            d_pair_pt_vec[patch_num].push_back({});
            tree.knnSearch(UPoint(patch, node), num_pts, idx_vec, distance_vec);
            for (const auto& idx_in_pts : idx_vec) d_pair_pt_vec[patch_num][i].push_back(pts[idx_in_pts]);
            ++i;
        }
    }
}
//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
