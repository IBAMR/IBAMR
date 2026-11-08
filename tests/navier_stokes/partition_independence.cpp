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

// Test that the velocity ghost values, the interpolation and spreading of the IB integrator, and the convective
// operator of the fluid solver do not depend on the patch layout or the number of processes. The program stops unless
// each printed value equals its reference value of the input to the tolerance of check_against_reference(), and each
// ghost value equals the value on the owning patch to 1e-12, relative to the larger of that value and 1.

#include <ibamr/ConvectiveOperator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LMesh.h>
#include <ibtk/LNode.h>
#include <ibtk/LNodeSetData.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
struct MarkerParameters
{
    int finest_ln;
    int num_cells;
    double region_upper;
};

// Place markers on lines parallel to each axis, a fraction of a grid cell inside the lower or upper boundary along each
// other axis, so that the kernels overlap physical boundaries, corners, and patch boundaries.
void
generate_markers(const unsigned int& /*strct_num*/,
                 const int& ln,
                 int& num_vertices,
                 std::vector<IBTK::Point>& vertex_posn,
                 void* ctx)
{
    const MarkerParameters& params = *static_cast<const MarkerParameters*>(ctx);
    vertex_posn.clear();
    if (ln == params.finest_ln)
    {
        const double h = 1.0 / params.num_cells;
        const int num_along_line = 2 * params.num_cells;
        int num_sides = 1;
        for (unsigned int d = 1; d < NDIM; ++d)
        {
            num_sides *= 2;
        }
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            // Each other axis contributes a lower or an upper offset.
            for (int sides = 0; sides < num_sides; ++sides)
            {
                for (int k = 0; k < num_along_line; ++k)
                {
                    IBTK::Point X;
                    int remaining_sides = sides;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            X[d] = (k + 0.5) / num_along_line;
                        }
                        else
                        {
                            const bool upper = remaining_sides % 2 == 1;
                            remaining_sides /= 2;
                            X[d] = upper ? 1.0 - (0.35 + 0.05 * d) * h : (0.3 + 0.05 * d) * h;
                        }
                    }
                    vertex_posn.push_back(X);
                }
            }
        }
        std::vector<IBTK::Point> kept;
        for (const IBTK::Point& X : vertex_posn)
        {
            bool inside = true;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                inside = inside && X[d] < params.region_upper;
            }
            if (inside)
            {
                kept.push_back(X);
            }
        }
        vertex_posn = kept;
    }
    num_vertices = static_cast<int>(vertex_posn.size());
    return;
} // generate_markers

// Move each marker by shift grid cells along each axis, without assigning the markers to patches again.
void
shift_markers(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
              Pointer<IBMethod> ib_method_ops,
              const std::array<double, NDIM>& shift)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = patch_hierarchy->getGridGeometry();
    const IntVector<NDIM>& ratio = patch_hierarchy->getPatchLevel(finest_ln)->getRatio();
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    Pointer<LData> X_data = l_data_manager->getLData("X", finest_ln);
    {
        boost::multi_array_ref<double, 2>& X = *X_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(finest_ln)->getLocalNodes())
        {
            const int petsc_idx = node->getLocalPETScIndex();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X[petsc_idx][d] += shift[d] * grid_geom->getDx()[d] / ratio(d);
            }
        }
        X_data->restoreArrays();
    }
    X_data->beginGhostUpdate();
    X_data->endGhostUpdate();
    return;
} // shift_markers

// The domain of the level, and the directions along which it is periodic.
struct LevelDomain
{
    Box<NDIM> box;
    std::array<bool, NDIM> periodic;
};

// Return whether the patch contains the cell c after c is mapped into the domain: to its periodic image along a
// periodic direction, and to the nearest cell of the domain along the others.
bool
patch_owns_cell(const Patch<NDIM>& patch, const LevelDomain& domain, hier::Index<NDIM> c)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const int lower = domain.box.lower(d), n = domain.box.numberCells(d);
        c(d) = domain.periodic[d] ? lower + ((c(d) - lower) % n + n) % n :
                                    std::max(lower, std::min(domain.box.upper(d), c(d)));
    }
    return patch.getBox().contains(c);
}

// Return whether the face i of the given axis lies beyond a physical boundary.
bool
is_outside(const LevelDomain& domain, const unsigned int axis, const hier::Index<NDIM>& i)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const int last = domain.box.upper(d) + (d == axis ? 1 : 0);
        if (!domain.periodic[d] && (i(d) < domain.box.lower(d) || i(d) > last)) return true;
    }
    return false;
}

// Return whether the face i lies in a periodic image of the domain, in the ghost cells of the patches, and not in the
// domain itself. The face at the upper end of a periodic direction is identified with the face at the lower end.
bool
is_periodic_image(const LevelDomain& domain, const hier::Index<NDIM>& i)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (domain.periodic[d] && (i(d) < domain.box.lower(d) || i(d) > domain.box.upper(d))) return true;
    }
    return false;
}

// Return the face of the domain that is the periodic image of the face i.
hier::Index<NDIM>
move_into_domain(const LevelDomain& domain, hier::Index<NDIM> i)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const int lower = domain.box.lower(d), n = domain.box.numberCells(d);
        if (domain.periodic[d]) i(d) = lower + ((i(d) - lower) % n + n) % n;
    }
    return i;
}

// Dense arrays over the faces of each axis within a ghost width of the domain.
class FaceArray
{
public:
    FaceArray(const LevelDomain& domain, const int ghost_width) : d_domain(domain), d_ghost_width(ghost_width)
    {
        d_size_per_axis = 1;
        for (unsigned int d = 0; d < NDIM; ++d) d_size_per_axis *= domain.box.numberCells(d) + 2 * ghost_width + 1;
        d_values.assign(NDIM * d_size_per_axis, 0.0);
        d_counts.assign(NDIM * d_size_per_axis, 0.0);
    }

    int flat_index(const unsigned int axis, const hier::Index<NDIM>& i) const
    {
        int k = 0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const int extent = d_domain.box.numberCells(d) + 2 * d_ghost_width + 1;
            k = k * extent + (i(d) - d_domain.box.lower(d) + d_ghost_width);
        }
        return static_cast<int>(axis * d_size_per_axis) + k;
    }

    void add(const unsigned int axis, const hier::Index<NDIM>& i, const double value)
    {
        d_values[flat_index(axis, i)] += value;
        d_counts[flat_index(axis, i)] += 1.0;
    }

    void sum_over_processes()
    {
        IBTK_MPI::sumReduction(d_values.data(), static_cast<int>(d_values.size()));
        IBTK_MPI::sumReduction(d_counts.data(), static_cast<int>(d_counts.size()));
    }

    std::vector<double> d_values, d_counts;

private:
    LevelDomain d_domain;
    int d_ghost_width;
    std::size_t d_size_per_axis;
};

// Print a quantity with enough digits to distinguish partitions, and return whether it equals the reference value of
// the input to 1e-12, relative to the largest of the reference value, 1, and term_size. The roundoff of a sum whose
// terms cancel scales with the root of the sum of the squares of its terms. The caller passes that root as term_size.
bool
check_against_reference(Pointer<Database> reference_db,
                        const std::string& key,
                        const std::string& label,
                        const double value,
                        const double term_size = 0.0)
{
    const double reference = reference_db->getDouble(key);
    const bool equal = std::abs(value - reference) <= 1.0e-12 * std::max({ 1.0, std::abs(reference), term_size });
    pout << label << " = " << value << '\n';
    if (!equal) pout << "  differs from the reference value " << reference << '\n';
    return equal;
} // check_against_reference

// Return the largest difference between a patch's copy of a value and the owning patch's copy, over the ghost locations
// beyond a physical boundary that the marker kernels of the patch use. The function stops with an error if a copy is
// not finite or differs from the owner's copy by more than 1e-12, relative to the larger of the owner's copy and 1.
double
largest_difference_from_owner_copy(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                   const LevelDomain& domain,
                                   const int u_idx,
                                   const int ghost_width,
                                   Pointer<LData> X_data,
                                   LDataManager* l_data_manager,
                                   const std::string& kernel)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);
    FaceArray owner_copies(domain, ghost_width);
    const double half_support = 0.5 * LEInteractor::getStencilSize(kernel);
    const double* const X = X_data->getGhostedLocalFormVecArray()->data();
    std::vector<int> local_indices;
    std::vector<double> periodic_shifts;
    struct Read
    {
        unsigned int axis;
        hier::Index<NDIM> i;
        double value;
    };
    std::vector<Read> reads;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
            {
                hier::Index<NDIM> c = b();
                c(axis) = std::min(c(axis), domain.box.upper(axis));
                if (is_outside(domain, axis, b()) && !is_periodic_image(domain, b()) &&
                    patch_owns_cell(*patch, domain, c))
                {
                    owner_copies.add(axis, b(), u_data->getArrayData(axis)(b(), 0));
                }
            }
        }

        // The faces in the kernel support of each marker in the patch box.
        Pointer<LNodeSetData> idx_data = patch->getPatchData(l_data_manager->getLNodePatchDescriptorIndex());
        LEInteractor::buildLocalIndices(local_indices, periodic_shifts, patch->getBox(), patch, X, NDIM, idx_data);
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        for (std::size_t n = 0; n < local_indices.size(); ++n)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                hier::Index<NDIM> lower, upper;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const double x = X[NDIM * local_indices[n] + d] + periodic_shifts[NDIM * n + d];
                    const double t = (x - pgeom->getXLower()[d]) / pgeom->getDx()[d] + patch->getBox().lower(d) -
                                     (d == axis ? 0.0 : 0.5);
                    lower(d) = static_cast<int>(std::floor(t - half_support)) + 1;
                    upper(d) = static_cast<int>(std::ceil(t + half_support)) - 1;
                }
                for (Box<NDIM>::Iterator b(Box<NDIM>(lower, upper)); b; b++)
                {
                    if (!is_outside(domain, axis, b())) continue;
                    if (!SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis).contains(b()))
                    {
                        TBOX_ERROR("the stencil of a marker includes the face "
                                   << b() << " of axis " << axis << ", beyond the ghost cells of its patch\n");
                    }
                    reads.push_back({ axis, b(), u_data->getArrayData(axis)(b(), 0) });
                }
            }
        }
    }
    X_data->restoreArrays();
    owner_copies.sum_over_processes();
    double largest = 0.0;
    for (const Read& read : reads)
    {
        const int flat_index = owner_copies.flat_index(read.axis, move_into_domain(domain, read.i));
        if (owner_copies.d_counts[flat_index] != 1.0)
        {
            TBOX_ERROR("the face " << read.i << " of axis " << read.axis << " has " << owner_copies.d_counts[flat_index]
                                   << " owner copies, not one\n");
        }
        const double difference = std::abs(read.value - owner_copies.d_values[flat_index]);
        if (!(difference <= 1.0e-12 * std::max(1.0, std::abs(owner_copies.d_values[flat_index]))))
        {
            TBOX_ERROR("the copy " << read.value << " of the face " << read.i << " of axis " << read.axis
                                   << " differs from the owner's copy " << owner_copies.d_values[flat_index] << '\n');
        }
        largest = std::max(largest, difference);
    }
    return IBTK_MPI::maxReduction(largest);
} // largest_difference_from_owner_copy

// Print sums and maxima of the interpolated velocity, the spread force density, the convective operator, and the
// velocity ghost values, counting each marker, face, and ghost face once, and compare them with the reference values of
// the input.
void
check_partition_independence(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                             Pointer<INSHierarchyIntegrator> ins_integrator,
                             Pointer<IBHierarchyIntegrator> ib_integrator,
                             Pointer<IBMethod> ib_method_ops,
                             const std::string& kernel,
                             Pointer<Database> reference_db)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);
    const double time = ib_integrator->getIntegratorTime();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const int u_idx = var_db->registerVariableAndContext(
        u_var, var_db->getContext("partition_independence"), ib_method_ops->getMinimumGhostCellWidth());
    const int f_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
    const int N_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        // The transfers use the data of every level that contains Lagrangian data.
        Pointer<PatchLevel<NDIM>> patch_level = patch_hierarchy->getPatchLevel(ln);
        patch_level->allocatePatchData(u_idx, time);
        patch_level->allocatePatchData(f_idx, time);
        patch_level->allocatePatchData(N_idx, time);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<SideData<NDIM, double>> u_data = patch_level->getPatch(p())->getPatchData(u_idx);
            Pointer<SideData<NDIM, double>> f_data = patch_level->getPatch(p())->getPatchData(f_idx);
            Pointer<SideData<NDIM, double>> N_data = patch_level->getPatch(p())->getPatchData(N_idx);
            u_data->fillAll(0.0);
            f_data->fillAll(0.0);
            N_data->fillAll(0.0);
        }
    }
    LevelDomain domain;
    BoxArray<NDIM> domain_boxes = patch_hierarchy->getGridGeometry()->getPhysicalDomain();
    domain_boxes.refine(level->getRatio());
    domain.box = domain_boxes[0];
    const IntVector<NDIM> periodic_shift = patch_hierarchy->getGridGeometry()->getPeriodicShift();
    for (unsigned int d = 0; d < NDIM; ++d) domain.periodic[d] = periodic_shift(d) != 0;

    // The velocity is a function of position, with the period of the domain along a periodic direction. The
    // divergence-free velocity boundary operator of the fluid solver sets the velocity ghost values outside the domain.
    RobinPhysBdryPatchStrategy* bdry_op = ib_integrator->getVelocityPhysBdryOp();
    bdry_op->setPatchDataIndex(u_idx);
    bdry_op->setHomogeneousBc(false);
    Pointer<RefineAlgorithm<NDIM>> ghost_fill_alg = new RefineAlgorithm<NDIM>();
    ghost_fill_alg->registerRefine(u_idx, u_idx, u_idx, nullptr);
    std::vector<Pointer<RefineSchedule<NDIM>>> ghost_fill_scheds(finest_ln + 1);
    ghost_fill_scheds[finest_ln] = ghost_fill_alg->createSchedule(level, bdry_op);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        u_data->fillAll(0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                double phase = 1.7 * axis;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const double x = x_lower[d] + dx[d] * (b()(d) - patch_box.lower(d) + (d == axis ? 0.0 : 0.5));
                    phase += (domain.periodic[d] ? 2.0 * M_PI * (1.0 + d) : 3.1 + 2.2 * d) * x;
                }
                u_data->getArrayData(axis)(b(), 0) = std::sin(phase);
            }
        }
    }
    ghost_fill_scheds[finest_ln]->fillData(time);

    // The velocity ghost values outside the domain, in the layers next to the patches that own them.
    double ghost_sum = 0.0, ghost_sum_squares = 0.0, ghost_max = 0.0, ghost_count = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
            {
                const hier::Index<NDIM>& i = b();
                hier::Index<NDIM> c = i;
                c(axis) = std::min(c(axis), domain.box.upper(axis));
                if (!is_outside(domain, axis, i) || is_periodic_image(domain, i) || !patch_owns_cell(*patch, domain, c))
                    continue;
                const double v = u_data->getArrayData(axis)(i, 0);
                ghost_sum += v;
                ghost_sum_squares += v * v;
                ghost_max = std::max(ghost_max, std::abs(v));
                ghost_count += 1.0;
            }
        }
    }

    // Interpolate u, and set a nodal force that depends only on the Lagrangian index.
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    std::vector<Pointer<LData>> X_vec(finest_ln + 1), U_vec(finest_ln + 1), F_vec(finest_ln + 1);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        X_vec[ln] = l_data_manager->getLData("X", ln);
        U_vec[ln] = l_data_manager->createLData("partition_independence_U", ln, NDIM);
        F_vec[ln] = l_data_manager->createLData("partition_independence_F", ln, NDIM);
    }
    l_data_manager->interp(u_idx, U_vec, X_vec, {}, ghost_fill_scheds, time, finest_ln, finest_ln);
    Pointer<LData> U_data = U_vec[finest_ln];
    Pointer<LData> F_data = F_vec[finest_ln];
    double F_dot_Ju = 0.0, F_dot_Ju_squares = 0.0, U_sum_squares = 0.0, U_max = 0.0, num_markers = 0.0;
    {
        boost::multi_array_ref<double, 2>& U = *U_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& F = *F_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(finest_ln)->getLocalNodes())
        {
            const int lag_idx = node->getLagrangianIndex();
            const int petsc_idx = node->getLocalPETScIndex();
            num_markers += 1.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F[petsc_idx][d] = std::cos(1.3 * lag_idx + 0.7 * d + 0.2);
                F_dot_Ju += F[petsc_idx][d] * U[petsc_idx][d];
                F_dot_Ju_squares += std::pow(F[petsc_idx][d] * U[petsc_idx][d], 2);
                U_sum_squares += U[petsc_idx][d] * U[petsc_idx][d];
                U_max = std::max(U_max, std::abs(U[petsc_idx][d]));
            }
        }
        U_data->restoreArrays();
        F_data->restoreArrays();
    }

    // Spread the nodal forces. The force density is zero outside the domain. Count each face once.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<SideData<NDIM, double>> f_data = level->getPatch(p())->getPatchData(f_idx);
        f_data->fillAll(0.0);
    }
    l_data_manager->spread(f_idx, F_vec, X_vec, bdry_op, {}, time, true, true, finest_ln, finest_ln);
    double f_sum = 0.0, f_sum_squares = 0.0, f_max = 0.0, SF_dot_u = 0.0, SF_dot_u_squares = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); b; b++)
            {
                // The face on the upper boundary of a periodic direction is the face on the lower boundary.
                if (domain.periodic[axis] && b()(axis) == domain.box.upper(axis) + 1) continue;
                hier::Index<NDIM> c = b();
                c(axis) = std::min(c(axis), domain.box.upper(axis));
                if (!patch_owns_cell(*patch, domain, c)) continue;
                const double v = f_data->getArrayData(axis)(b(), 0);
                f_sum += v;
                f_sum_squares += v * v;
                f_max = std::max(f_max, std::abs(v));
                SF_dot_u += v * u_data->getArrayData(axis)(b(), 0);
                SF_dot_u_squares += std::pow(v * u_data->getArrayData(axis)(b(), 0), 2);
            }
        }
    }

    // Apply the convective operator of the fluid solver to u.
    Pointer<ConvectiveOperator> N_op = ins_integrator->getConvectiveOperator();
    N_op->setAdvectionVelocity(u_idx);
    N_op->setSolutionTime(time);
    SAMRAIVectorReal<NDIM, double> u_vec("partition_independence_u", patch_hierarchy, 0, finest_ln);
    SAMRAIVectorReal<NDIM, double> N_vec("partition_independence_N", patch_hierarchy, 0, finest_ln);
    u_vec.addComponent(u_var, u_idx, -1);
    N_vec.addComponent(u_var, N_idx, -1);
    N_op->initializeOperatorState(u_vec, N_vec);
    N_op->applyConvectiveOperator(u_idx, N_idx);
    N_op->deallocateOperatorState();
    double N_sum = 0.0, N_sum_squares = 0.0, N_max = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> N_data = patch->getPatchData(N_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); b; b++)
            {
                if (domain.periodic[axis] && b()(axis) == domain.box.upper(axis) + 1) continue;
                hier::Index<NDIM> c = b();
                c(axis) = std::min(c(axis), domain.box.upper(axis));
                if (!patch_owns_cell(*patch, domain, c)) continue;
                const double v = N_data->getArrayData(axis)(b(), 0);
                N_sum += v;
                N_sum_squares += v * v;
                N_max = std::max(N_max, std::abs(v));
            }
        }
    }

    const double copy_difference = largest_difference_from_owner_copy(patch_hierarchy,
                                                                      domain,
                                                                      u_idx,
                                                                      ib_method_ops->getMinimumGhostCellWidth().max(),
                                                                      X_vec[finest_ln],
                                                                      l_data_manager,
                                                                      kernel);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(u_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(f_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(N_idx);
    }

    pout << std::setprecision(14);
    bool all_equal = true;
    all_equal &= check_against_reference(reference_db, "markers", "markers", IBTK_MPI::sumReduction(num_markers));
    all_equal &=
        check_against_reference(reference_db, "ghost_faces", "ghost faces", IBTK_MPI::sumReduction(ghost_count));
    all_equal &= check_against_reference(reference_db,
                                         "ghost_sum",
                                         "sum of ghost values",
                                         IBTK_MPI::sumReduction(ghost_sum),
                                         std::sqrt(IBTK_MPI::sumReduction(ghost_sum_squares)));
    all_equal &= check_against_reference(
        reference_db, "ghost_sum_squares", "sum of squares of ghost values", IBTK_MPI::sumReduction(ghost_sum_squares));
    all_equal &=
        check_against_reference(reference_db, "ghost_max", "max ghost value", IBTK_MPI::maxReduction(ghost_max));
    pout << "largest difference between a ghost copy and the owner's copy = " << copy_difference << '\n';
    all_equal &= check_against_reference(reference_db,
                                         "F_dot_Ju",
                                         "(F, J u)",
                                         IBTK_MPI::sumReduction(F_dot_Ju),
                                         std::sqrt(IBTK_MPI::sumReduction(F_dot_Ju_squares)));
    all_equal &= check_against_reference(
        reference_db, "Ju_sum_squares", "sum of squares of J u", IBTK_MPI::sumReduction(U_sum_squares));
    all_equal &= check_against_reference(reference_db, "Ju_max", "max J u", IBTK_MPI::maxReduction(U_max));
    all_equal &= check_against_reference(reference_db,
                                         "SF_sum",
                                         "sum of S F",
                                         IBTK_MPI::sumReduction(f_sum),
                                         std::sqrt(IBTK_MPI::sumReduction(f_sum_squares)));
    all_equal &= check_against_reference(
        reference_db, "SF_sum_squares", "sum of squares of S F", IBTK_MPI::sumReduction(f_sum_squares));
    all_equal &= check_against_reference(reference_db, "SF_max", "max S F", IBTK_MPI::maxReduction(f_max));
    all_equal &= check_against_reference(reference_db,
                                         "SF_dot_u",
                                         "(S F, u)",
                                         IBTK_MPI::sumReduction(SF_dot_u),
                                         std::sqrt(IBTK_MPI::sumReduction(SF_dot_u_squares)));
    all_equal &= check_against_reference(reference_db,
                                         "N_sum",
                                         "sum of N",
                                         IBTK_MPI::sumReduction(N_sum),
                                         std::sqrt(IBTK_MPI::sumReduction(N_sum_squares)));
    all_equal &= check_against_reference(
        reference_db, "N_sum_squares", "sum of squares of N", IBTK_MPI::sumReduction(N_sum_squares));
    all_equal &= check_against_reference(reference_db, "N_max", "max N", IBTK_MPI::maxReduction(N_max));
    if (!all_equal) TBOX_ERROR("the printed quantities differ from the reference values\n");
    return;
} // check_partition_independence
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Suppress warnings caused by running without Silo.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        Pointer<IBHierarchyIntegrator> ib_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              ins_integrator);
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ib_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        MarkerParameters marker_params{ input_db->getInteger("MAX_LEVELS") - 1,
                                        input_db->getInteger("N"),
                                        input_db->getDoubleWithDefault("marker_region_upper", 1.0) };
        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        ib_initializer->setStructureNamesOnLevel(marker_params.finest_ln, { "markers" });
        ib_initializer->registerInitStructureFunction(generate_markers, &marker_params);
        ib_method_ops->registerLInitStrategy(ib_initializer);
        ib_method_ops->registerIBLagrangianForceFunction(new IBStandardForceGen());

        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        ins_integrator->registerVelocityInitialConditions(u_init);
        std::vector<std::unique_ptr<muParserRobinBcCoefs>> u_bc_coef_owners;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_coef_owners.push_back(std::make_unique<muParserRobinBcCoefs>(
                "u_bc_coefs_" + std::to_string(d),
                app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                grid_geometry));
            u_bc_coefs.push_back(u_bc_coef_owners.back().get());
        }
        ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        ib_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();

        // Regrid, as advanceHierarchy() does at the initial time, to set up the Lagrangian indices in ghost cells.
        // Preprocessing sets up the velocity boundary conditions.
        ib_integrator->regridHierarchy();
        const double current_time = ib_integrator->getIntegratorTime();
        ib_integrator->preprocessIntegrateHierarchy(
            current_time, current_time + ib_integrator->getMaximumTimeStepSize(), 1);
        if (input_db->keyExists("marker_shift"))
        {
            std::array<double, NDIM> marker_shift;
            input_db->getDoubleArray("marker_shift", marker_shift.data(), NDIM);
            shift_markers(patch_hierarchy, ib_method_ops, marker_shift);
        }
        check_partition_independence(patch_hierarchy,
                                     ins_integrator,
                                     ib_integrator,
                                     ib_method_ops,
                                     app_initializer->getComponentDatabase("IBMethod")->getString("delta_fcn"),
                                     input_db->getDatabase("Reference"));
    }
} // main
