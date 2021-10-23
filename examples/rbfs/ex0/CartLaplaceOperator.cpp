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
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"

#include "CellVariable.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "tbox/Timer.h"

#include <Eigen/Dense>

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

static unsigned int interp_size = 22; // 2 * (NDIM + 1) + 1;

} // namespace

unsigned int CartLaplaceOperator::s_num_ghost_cells = 3;
/////////////////////////////// PUBLIC ///////////////////////////////////////

CartLaplaceOperator::CartLaplaceOperator(const std::string& object_name,
                                         BoundaryMesh* bdry_mesh,
                                         const DofMap* dof_map,
                                         Pointer<Database> input_db)
    : PETScLinearAugmentedOperator(object_name, false), d_mesh(bdry_mesh), d_dof_map(dof_map)
{
    d_dist_to_bdry = input_db->getDouble("dist_to_bdry");
    d_eps = input_db->getDouble("eps");
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
#if (0)
        if (x_depth != d_bc_coefs.size() || y_depth != d_bc_coefs.size())
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  each vector component must have data depth == " << d_bc_coefs.size() << "\n"
                                     << "  since d_bc_coefs.size() == " << d_bc_coefs.size() << std::endl);
        }
#endif
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
                                                      CONSISTENT_TYPE_2_BDRY);
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

    applyToLagDOFs(x.getComponentDescriptorIndex(0), y.getComponentDescriptorIndex(0));

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
                                                    CONSISTENT_TYPE_2_BDRY);
        d_transaction_comps.push_back(component);
    }

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy, d_coarsest_ln, d_finest_ln);

    // Sort Lag DOFs
    sortLagDOFsToCells();

    // Clone Lag Vector
    int ierr = VecDuplicate(d_aug_x_vec, &d_aug_b_vec);
    IBTK_CHKERRQ(ierr);

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
    //    auto rbf = [](const double r) -> double { return r * r * r; };
    //    auto lap_rbf = [](const double r) -> double {return 9.0 * r;};
    auto rbf = [](const double r) -> double { return r * r * r * r * r + 2.0e-10; };
    auto lap_rbf = [](const double r) -> double { return 25.0 * r * r * r; };
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_num)
    {
        if (d_base_pt_vec[patch_num].size() == 0) continue;
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        // First loop through Cartesian grid cells.
        const Box<NDIM>& box = patch->getBox();
        Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(d_ls_idx);
        Pointer<CellData<NDIM, double> > x_data = patch->getPatchData(x_idx);
        Pointer<CellData<NDIM, double> > y_data = patch->getPatchData(y_idx);
        // Note Lagrangian data are located in d_aug_(x|b)_vec

        // All data have been sorted. We need to loop through d_base_pt_vec.
        for (size_t idx = 0; idx < d_base_pt_vec[patch_num].size(); ++idx)
        {
            const UPoint& pt = d_base_pt_vec[patch_num][idx];
            const std::vector<UPoint>& pt_vec = d_pair_pt_vec[patch_num][idx];
            // Note if we use a KNN search, interp_size is fixed.
            const int interp_size = pt_vec.size();
            // Up to cubic polynomials
            const int poly_size = NDIM + 1 + NDIM + 1 + NDIM * NDIM;
#define DEBUGGING 1
#if (DEBUGGING)
            plog << "On point \n" << pt << "\n";
            plog << "Forming interpolant with " << interp_size << " points and " << poly_size << " polynomials\n";
#endif
            MatrixXd A(MatrixXd::Zero(interp_size, interp_size));
            MatrixXd B(MatrixXd::Zero(interp_size, poly_size));
            VectorXd U(VectorXd::Zero(interp_size + poly_size));
            VectorNd pt0 = pt.getVec();
            for (int d = 0; d < NDIM; ++d) pt0[d] = pt0[d] / dx[d];
            for (size_t i = 0; i < interp_size; ++i)
            {
                VectorNd pti = pt_vec[i].getVec();
                for (int d = 0; d < NDIM; ++d) pti[d] = pti[d] / dx[d];
                for (size_t j = 0; j < interp_size; ++j)
                {
                    VectorNd ptj = pt_vec[j].getVec();
                    for (int d = 0; d < NDIM; ++d) ptj[d] = ptj[d] / dx[d];
                    A(i, j) = rbf((pti - ptj).norm());
                }
                // TODO: B is just a Vandermonde matrix. Write a function to set this up given arbitrary polynomial
                // degree.
                B(i, 0) = 1.0;
                VectorNd diff = pti - pt0;
                for (int d = 0; d < NDIM; ++d) B(i, d + 1) = diff(d);
                // Add quadratic polynomials
                B(i, NDIM + 1) = diff(0) * diff(0);
                B(i, NDIM + 2) = diff(1) * diff(0);
                B(i, NDIM + 3) = diff(1) * diff(1);
                // Cubic
                B(i, NDIM + 4) = diff(0) * diff(0) * diff(0);
                B(i, NDIM + 5) = diff(1) * diff(0) * diff(0);
                B(i, NDIM + 6) = diff(1) * diff(1) * diff(0);
                B(i, NDIM + 7) = diff(1) * diff(1) * diff(1);
                // Determine rhs
                U(i) = lap_rbf((pt0 - pti).norm());
            }
            // Add quadratic polynomials
            U(interp_size + NDIM + 1) = 2.0;
            U(interp_size + NDIM + 3) = 2.0;
            MatrixXd final_mat(MatrixXd::Zero(interp_size + poly_size, interp_size + poly_size));
            final_mat.block(0, 0, interp_size, interp_size) = A;
            final_mat.block(0, interp_size, interp_size, poly_size) = B;
            final_mat.block(interp_size, 0, poly_size, interp_size) = B.transpose();

            VectorXd x = final_mat.fullPivHouseholderQr().solve(U);
#if (DEBUGGING)
            plog << "After assembly, A is:\n " << A << "\n";
            plog << "B is: \n" << B << "\n";
            plog << "U is : \n" << U << "\n";
            plog << "final mat: \n" << final_mat << "\n";
            plog << "final mat det: " << final_mat.determinant() << "\n";
            plog << "Residual: \n" << final_mat * x - U << "\n";
            plog << "Solution is : \n" << x << "\n";
#endif
            // Now evaluate FD stencil)
            double val = 0.0;
            VectorXd weights = x.block(0, 0, interp_size, 1);
            for (size_t i = 0; i < interp_size; ++i)
            {
                double w = weights[i];
#if (DEBUGGING)
                plog << "Solution value at point " << i << " is " << getSolVal(pt_vec[i], *x_data) << "\n";
                plog << "Weight at point         " << i << " is " << w << "\n";
#endif
                val += w * getSolVal(pt_vec[i], *x_data) / (dx[0] * dx[1]);
            }
            // Now insert val into results
            setSolVal(val, pt, *y_data);
#if (DEBUGGING)
            plog << "\n";
#endif
        }
    }
}

void
CartLaplaceOperator::sortLagDOFsToCells()
{
    // Fill ghost cells for level set.
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_comps;
    ghost_cell_comps.push_back(ITC(
        d_ls_idx, DATA_REFINE_TYPE, USE_CF_INTERPOLATION, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY));
    HierarchyGhostCellInterpolation ghost_cell_fill;
    ghost_cell_fill.initializeOperatorState(ghost_cell_comps, d_hierarchy);
    ghost_cell_fill.fillData(0.0);

    // Clear old data structures.
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
        // We also need ghost nodes to include in our KD tree
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
        // We are on a patch that has points. We need to form a KD tree.
        // TODO: Need to determine when we need to use RBFs
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(d_ls_idx);

        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        // First start by collecting all points into a vector
        std::vector<UPoint> pts;
        for (CellIterator<NDIM> ci(ls_data->getGhostBox()); ci; ci++)
        {
            const CellIndex<NDIM>& idx = ci();
            if ((*ls_data)(idx) < -d_eps) pts.push_back(UPoint(patch, idx));
        }
        // Now add in the points in d_idx_node_vec
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
            if (ls < -d_eps && std::abs(ls) < d_dist_to_bdry)
            {
                std::vector<int> idx_vec;
                std::vector<double> distance_vec;
                d_base_pt_vec[patch_num].push_back(UPoint(patch, idx));
                d_pair_pt_vec[patch_num].push_back({});
#if (1)
                // Use KNN search.
                tree.knnSearch(UPoint(patch, idx), interp_size, idx_vec, distance_vec);
#else
                // Use a bounding box search
                VectorNd bbox;
                bbox(0) = 2.0 * dx[0];
                bbox(1) = 2.0 * dx[1];
                tree.cuboid_query(UPoint(patch, idx), bbox, idx_vec, distance_vec);
#endif
                // Add these points to the vector
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
            tree.knnSearch(UPoint(patch, node), interp_size, idx_vec, distance_vec);
            // Add these points to the vector
            for (const auto& idx_in_pts : idx_vec) d_pair_pt_vec[patch_num][i].push_back(pts[idx_in_pts]);
            ++i;
        }
    }
}

double
CartLaplaceOperator::getSolVal(const UPoint& pt, const CellData<NDIM, double>& Q_data) const
{
    double val = 0.0;
    if (pt.isNode())
    {
        // We're on a node. Need to grab value from augmented vec
        std::vector<unsigned int> dof_indices;
        d_dof_map->dof_indices(pt.getNode(), dof_indices);
        auto idxs = reinterpret_cast<PetscInt*>(dof_indices.data());
        int ierr = VecGetValues(d_aug_x_vec, 1, idxs, &val);
        IBTK_CHKERRQ(ierr);
    }
#ifndef NDEBUG
    else if (pt.isEmpty())
    {
        TBOX_ERROR("Point should not be empty");
    }
#endif
    else
    {
        val = Q_data(pt.getIndex());
    }
    return val;
}

void
CartLaplaceOperator::setSolVal(const double val, const UPoint& pt, CellData<NDIM, double>& Q_data) const
{
    if (pt.isNode())
    {
        // We're on a node. Need to grab value from augmented vec
        std::vector<unsigned int> dof_indices;
        d_dof_map->dof_indices(pt.getNode(), dof_indices);
        int ierr = VecSetValue(d_aug_b_vec, dof_indices[0], val, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
#ifndef NDEBUG
    else if (pt.isEmpty())
    {
        TBOX_ERROR("Point should not be empty");
    }
#endif
    else
    {
        Q_data(pt.getIndex()) = val;
    }
}
//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
