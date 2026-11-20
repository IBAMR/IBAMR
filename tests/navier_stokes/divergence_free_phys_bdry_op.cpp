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

#include <ibamr/INSHierarchyIntegrator.h>
#include <ibamr/INSStaggeredDivergenceFreePhysBdryOp.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>
#include <muParser.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{

// A velocity given by muParser expressions function_0, function_1, ... of the
// position X_0, X_1, ....
class ParsedVelocity
{
public:
    explicit ParsedVelocity(Pointer<Database> db) : d_parsers(NDIM)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_parsers[d].SetExpr(db->getString("function_" + std::to_string(d)));
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                d_parsers[d].DefineVar("X_" + std::to_string(k), &d_X[k]);
            }
        }
    }

    ParsedVelocity(const ParsedVelocity&) = delete;
    ParsedVelocity& operator=(const ParsedVelocity&) = delete;

    double operator()(const unsigned int axis, const std::array<double, NDIM>& x)
    {
        d_X = x;
        return d_parsers[axis].Eval();
    }

private:
    std::array<double, NDIM> d_X;
    std::vector<mu::Parser> d_parsers;
};

void
check_divergence_free_extension(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                Pointer<INSHierarchyIntegrator> ins_integrator,
                                Pointer<Database> exact_velocity_db,
                                const int extension_ghost_width,
                                const int u_idx,
                                const int y_idx)
{
    // Fill a multi-patch level and check the ghost divergence, the transpose identity (y, E u) = (E^T y, u) of the
    // extension E, and the error in reproducing the velocity.
    if (IBTK_MPI::getNodes() > 1 || patch_hierarchy->getFinestLevelNumber() > 0)
    {
        TBOX_ERROR("check_divergence_free_extension():\n"
                   << "  this check requires one rank and one level, not " << IBTK_MPI::getNodes() << " ranks and "
                   << patch_hierarchy->getFinestLevelNumber() + 1 << " levels\n"
                   << "  run on one rank and set max_levels = 1 in the input\n");
    }
    const int G = extension_ghost_width;
    ParsedVelocity exact_velocity(exact_velocity_db);
    const IntVector<NDIM> ghost_width(G);
    INSStaggeredDivergenceFreePhysBdryOp bdry_op(ins_integrator, /*homogeneous_bc*/ false);
    INSStaggeredDivergenceFreePhysBdryOp homogeneous_bdry_op(ins_integrator, /*homogeneous_bc*/ true);
    const double fill_time = ins_integrator->getIntegratorTime();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const bool has_source = ins_integrator->hasVelocityDivergenceFunction();
    const bool approximate = has_source || exact_velocity_db->keyExists("max_error_over_h_squared");
    double max_dx = 0.0;
    {
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = patch_hierarchy->getGridGeometry();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            max_dx = std::max(max_dx, grid_geometry->getDx()[d]);
        }
    }

    // A deterministic function of the face that is not smooth.
    const auto rough = [](const unsigned int axis, const hier::Index<NDIM>& i, const double scale)
    {
        static const std::array<double, 3> coefficients = { 12.9898, 78.233, 37.719 };
        double phase = 4.1 * axis;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            phase += coefficients[d] * i(d);
        }
        return std::sin(scale * phase);
    };

    double max_divergence = 0.0;
    double divergence_scale = 0.0;
    double patch_difference = 0.0;
    double ghost_scale = 0.0;
    double y_dot_Eu = 0.0;
    double ETy_dot_u = 0.0;
    double max_filled_after_transpose = 0.0;
    double num_nonfinite_divergence = 0.0;
    double exact_error = 0.0;
    double num_nan = 0.0;
    std::map<std::array<int, NDIM + 1>, double> face_values;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_idx, fill_time);
        level->allocatePatchData(y_idx, fill_time);
        BoxArray<NDIM> domain = patch_hierarchy->getGridGeometry()->getPhysicalDomain();
        domain.refine(level->getRatio());
        const Box<NDIM> domain_box = domain[0];
        const IntVector<NDIM> periodic_shift = patch_hierarchy->getGridGeometry()->getPeriodicShift(level->getRatio());

        // Whether the face i of component axis (the cell i if axis = NDIM) lies beyond a boundary of the domain.
        const auto is_beyond_boundary = [&](const unsigned int axis, const hier::Index<NDIM>& i)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (periodic_shift(d) == 0 &&
                    (i(d) < domain_box.lower(d) || i(d) > domain_box.upper(d) + (d == axis ? 1 : 0)))
                {
                    return true;
                }
            }
            return false;
        };

        Pointer<RefineAlgorithm<NDIM>> ghost_fill_alg = new RefineAlgorithm<NDIM>();
        ghost_fill_alg->registerRefine(u_idx, u_idx, u_idx, nullptr);
        Pointer<RefineSchedule<NDIM>> ghost_fill_sched = ghost_fill_alg->createSchedule(level, &bdry_op);
        Pointer<RefineSchedule<NDIM>> homogeneous_ghost_fill_sched =
            ghost_fill_alg->createSchedule(level, &homogeneous_bdry_op);

        // The distance of the face i of component axis from the patch along each direction parallel to every boundary
        // beyond the face. The values within G - 1 of the patch are independent of the patch.
        const auto parallel_offset = [&](const Patch<NDIM>& patch, const unsigned int axis, const hier::Index<NDIM>& i)
        {
            int offset = 0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int last = domain_box.upper(d) + (d == axis ? 1 : 0);
                if (periodic_shift(d) == 0 && (i(d) < domain_box.lower(d) || i(d) > last))
                {
                    continue;
                }
                offset = std::max(
                    offset,
                    std::max(patch.getBox().lower(d) - i(d), i(d) - patch.getBox().upper(d) - (d == axis ? 1 : 0)));
            }
            return offset;
        };

        // Call fcn(patch, axis, i) for each side index i within width of the patch that is outside the domain.
        const auto for_each_exterior_face = [&](Pointer<Patch<NDIM>> patch, const int width, const auto& fcn)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> box = SideGeometry<NDIM>::toSideBox(Box<NDIM>::grow(patch->getBox(), width), axis);
                for (Box<NDIM>::Iterator b(box); b; b++)
                {
                    if (is_beyond_boundary(axis, b()))
                    {
                        fcn(axis, b());
                    }
                }
            }
        };

        // Set the values in the patch interiors and NaN elsewhere, so that a value that depends on an unset value is
        // NaN.
        const auto fill_ghosts = [&](const auto& value, const Pointer<RefineSchedule<NDIM>>& sched)
        {
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
                u_data->fillAll(nan);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); b; b++)
                    {
                        u_data->getArrayData(axis)(b(), 0) = value(*patch, axis, b());
                    }
                }
            }
            sched->fillData(fill_time);
        };

        // Update the maximum error of the ghost divergence, whose expected value is zero or, with a source, the mean of
        // the divergence of the cells shifted into the domain, and the scale max |u| / dx.
        const auto update_exterior_divergence = [&]()
        {
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
                for (Box<NDIM>::Iterator b(Box<NDIM>::grow(patch->getBox(), G)); b; b++)
                {
                    if (!is_beyond_boundary(NDIM, b()))
                    {
                        continue;
                    }
                    bool has_nan_face = false;
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        hier::Index<NDIM> i_upper = b();
                        i_upper(axis) += 1;
                        has_nan_face = has_nan_face || std::isnan(u_data->getArrayData(axis)(b(), 0)) ||
                                       std::isnan(u_data->getArrayData(axis)(i_upper, 0));
                    }
                    if (has_nan_face)
                    {
                        continue;
                    }
                    const auto cell_divergence = [&](const hier::Index<NDIM>& cell)
                    {
                        double sum = 0.0;
                        for (unsigned int axis = 0; axis < NDIM; ++axis)
                        {
                            hier::Index<NDIM> i_upper = cell;
                            i_upper(axis) += 1;
                            const double u_upper = u_data->getArrayData(axis)(i_upper, 0);
                            const double u_lower = u_data->getArrayData(axis)(cell, 0);
                            sum += (u_upper - u_lower) / pgeom->getDx()[axis];
                            divergence_scale =
                                std::max(divergence_scale,
                                         std::max(std::abs(u_upper), std::abs(u_lower)) / pgeom->getDx()[axis]);
                        }
                        return sum;
                    };
                    double divergence = cell_divergence(b());
                    if (has_source)
                    {
                        double expected = 0.0;
                        int num_boundaries = 0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            if (periodic_shift(d) != 0 ||
                                (b()(d) >= domain_box.lower(d) && b()(d) <= domain_box.upper(d)))
                            {
                                continue;
                            }
                            hier::Index<NDIM> nearest = b();
                            nearest(d) = b()(d) < domain_box.lower(d) ? domain_box.lower(d) : domain_box.upper(d);
                            expected += cell_divergence(nearest);
                            ++num_boundaries;
                        }
                        divergence -= expected / num_boundaries;
                    }
                    if (std::isfinite(divergence))
                    {
                        max_divergence = std::max(max_divergence, std::abs(divergence));
                    }
                    else
                    {
                        num_nonfinite_divergence += 1.0;
                    }
                }
            }
        };

        const auto rough_value = [&](const Patch<NDIM>&, const unsigned int axis, const hier::Index<NDIM>& i)
        { return rough(axis, i, 1.0); };

        // Linear divergence-free velocity and boundary data.
        const auto exact_value = [&](const Patch<NDIM>& patch, const unsigned int axis, const hier::Index<NDIM>& i)
        {
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
            std::array<double, NDIM> X;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X[d] = pgeom->getXLower()[d] +
                       pgeom->getDx()[d] * (i(d) - patch.getBox().lower(d) + (d == axis ? 0.0 : 0.5));
            }
            return exact_velocity(axis, X);
        };
        bdry_op.setPatchDataIndex(u_idx);
        fill_ghosts(exact_value, ghost_fill_sched);
        update_exterior_divergence();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            for_each_exterior_face(patch,
                                   G,
                                   [&](const unsigned int axis, const hier::Index<NDIM>& i)
                                   {
                                       const double value = u_data->getArrayData(axis)(i, 0);
                                       if (std::isnan(value))
                                       {
                                           if (parallel_offset(*patch, axis, i) < G)
                                           {
                                               num_nan += 1.0;
                                           }
                                           return;
                                       }
                                       std::array<int, NDIM + 1> key;
                                       key[0] = static_cast<int>(axis);
                                       for (unsigned int d = 0; d < NDIM; ++d)
                                       {
                                           key[d + 1] = i(d);
                                       }
                                       ghost_scale = std::max(ghost_scale, std::abs(value));
                                       if (parallel_offset(*patch, axis, i) < G)
                                       {
                                           const auto it = face_values.find(key);
                                           if (it == face_values.end())
                                           {
                                               face_values[key] = value;
                                           }
                                           else
                                           {
                                               patch_difference =
                                                   std::max(patch_difference, std::abs(value - it->second));
                                           }
                                       }
                                       exact_error =
                                           std::max(exact_error, std::abs(value - exact_value(*patch, axis, i)));
                                   });
        }

        // Velocity that is not smooth.
        fill_ghosts(rough_value, ghost_fill_sched);
        update_exterior_divergence();

        // Transpose identity, with arbitrary y outside the domain within width G and u in the domain not smooth.
        homogeneous_bdry_op.setPatchDataIndex(u_idx);
        fill_ghosts(rough_value, homogeneous_ghost_fill_sched);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<SideData<NDIM, double>> y_data = patch->getPatchData(y_idx);
            y_data->fillAll(0.0);
            for_each_exterior_face(patch,
                                   G,
                                   [&](const unsigned int axis, const hier::Index<NDIM>& i)
                                   {
                                       if (std::isnan(u_data->getArrayData(axis)(i, 0))) return;
                                       y_data->getArrayData(axis)(i, 0) = rough(axis, i, 1.3) + 0.2;
                                       y_dot_Eu += y_data->getArrayData(axis)(i, 0) * u_data->getArrayData(axis)(i, 0);
                                   });
            homogeneous_bdry_op.setPatchDataIndex(y_idx);
            homogeneous_bdry_op.accumulateFromPhysicalBoundaryData(*patch, fill_time, ghost_width);
            homogeneous_bdry_op.setPatchDataIndex(u_idx);

            // The accumulation sets to zero every computable value outside the domain within the ghost width.
            for_each_exterior_face(patch,
                                   G,
                                   [&](const unsigned int axis, const hier::Index<NDIM>& i)
                                   {
                                       const double value = y_data->getArrayData(axis)(i, 0);
                                       if (!std::isfinite(value))
                                       {
                                           TBOX_ERROR("check_divergence_free_extension():\n"
                                                      << "  the value after the accumulation of component " << axis
                                                      << " at the face " << i << " outside the domain is not finite\n");
                                       }
                                       max_filled_after_transpose =
                                           std::max(max_filled_after_transpose, std::abs(value));
                                   });
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    if (!is_beyond_boundary(axis, b()))
                    {
                        ETy_dot_u += y_data->getArrayData(axis)(b(), 0) * u_data->getArrayData(axis)(b(), 0);
                    }
                }
            }
        }
        level->deallocatePatchData(u_idx);
        level->deallocatePatchData(y_idx);
    }

    const double divergence = IBTK_MPI::maxReduction(max_divergence);
    const double divergence_bound = 1.0e-10 * IBTK_MPI::maxReduction(divergence_scale);
    if (!(divergence <= divergence_bound))
    {
        TBOX_ERROR("check_divergence_free_extension():\n"
                   << "  the divergence in a ghost cell outside the domain is " << divergence << ", which exceeds "
                   << divergence_bound << "\n");
    }
    if (IBTK_MPI::sumReduction(num_nonfinite_divergence) != 0.0 || IBTK_MPI::sumReduction(num_nan) != 0.0)
    {
        TBOX_ERROR("check_divergence_free_extension():\n"
                   << "  a ghost value within width " << G << " is not finite\n");
    }
    const double patch_difference_max = IBTK_MPI::maxReduction(patch_difference);
    if (!(patch_difference_max <= 1.0e-12 * IBTK_MPI::maxReduction(ghost_scale)))
    {
        TBOX_ERROR("check_divergence_free_extension():\n"
                   << "  two patches differ by " << patch_difference_max << " in a ghost value\n");
    }
    if (IBTK_MPI::maxReduction(max_filled_after_transpose) != 0.0)
    {
        TBOX_ERROR("check_divergence_free_extension():\n"
                   << "  a value outside the domain is not zero after the accumulation\n");
    }
    const double y_dot_Eu_sum = IBTK_MPI::sumReduction(y_dot_Eu);
    const double ETy_dot_u_sum = IBTK_MPI::sumReduction(ETy_dot_u);
    if (!IBTK::rel_equal_eps(y_dot_Eu_sum, ETy_dot_u_sum, 1.0e-12))
    {
        TBOX_ERROR("check_divergence_free_extension():\n"
                   << "  (y, E u) = " << y_dot_Eu_sum << " and (E^T y, u) = " << ETy_dot_u_sum << " differ\n");
    }
    const double exact_error_max = IBTK_MPI::maxReduction(exact_error);
    pout << std::setprecision(12) << "(y, E u)   = " << y_dot_Eu_sum << '\n'
         << "(E^T y, u) = " << ETy_dot_u_sum << '\n';
    if (approximate)
    {
        // The error is a multiple of h^2 with a coefficient independent of h.
        const double error_over_h2 = exact_error_max / (max_dx * max_dx);
        if (!(error_over_h2 <= exact_velocity_db->getDouble("max_error_over_h_squared")))
        {
            TBOX_ERROR("check_divergence_free_extension():\n"
                       << "  the extension differs from the exact velocity by " << exact_error_max << ", which is "
                       << error_over_h2 << " h^2\n"
                       << "  the error over h^2 must not exceed max_error_over_h_squared in the input\n");
        }
        pout << std::setprecision(6) << "largest ghost value error / h^2 = " << error_over_h2 << '\n';
    }
    else if (!(exact_error_max <= 1.0e-12 * IBTK_MPI::maxReduction(ghost_scale)))
    {
        TBOX_ERROR("check_divergence_free_extension():\n"
                   << "  the extension of a linear divergence-free velocity differs from the exact velocity by "
                   << exact_error_max << "\n");
    }
    return;
} // check_divergence_free_extension

// Check that a fill of ghost width w less than the ghost width G of the data is independent of the unfilled faces
// inside the domain: two different finite values in them give identical results in every face of the data.
void
check_unfilled_data_unread(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                           Pointer<INSHierarchyIntegrator> ins_integrator,
                           const int fill_ghost_width,
                           const int u_idx)
{
    if (IBTK_MPI::getNodes() > 1 || patch_hierarchy->getFinestLevelNumber() > 0)
    {
        TBOX_ERROR("check_unfilled_data_unread():\n"
                   << "  this check requires one rank and one level\n"
                   << "  run on one rank and set max_levels = 1 in the input\n");
    }
    INSStaggeredDivergenceFreePhysBdryOp bdry_op(ins_integrator, /*homogeneous_bc*/ false);
    const double fill_time = ins_integrator->getIntegratorTime();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    level->allocatePatchData(u_idx, fill_time);
    const Box<NDIM> domain_box = patch_hierarchy->getGridGeometry()->getPhysicalDomain()[0];
    const IntVector<NDIM> periodic_shift = patch_hierarchy->getGridGeometry()->getPeriodicShift(level->getRatio());
    const auto is_beyond_boundary = [&](const unsigned int axis, const hier::Index<NDIM>& i)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (periodic_shift(d) == 0 &&
                (i(d) < domain_box.lower(d) || i(d) > domain_box.upper(d) + (d == axis ? 1 : 0)))
            {
                return true;
            }
        }
        return false;
    };
    bdry_op.setPatchDataIndex(u_idx);
    int num_different = 0, num_finite_outside = 0, num_nan_in_fill_box = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        std::array<std::vector<double>, 2> results;
        for (int run = 0; run < 2; ++run)
        {
            // Both runs set the same values in the patch and in the faces inside the domain within the fill width.
            // The unfilled faces inside the domain contain a value that differs between the runs.
            u_data->fillAll(nan);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> fill_box =
                    SideGeometry<NDIM>::toSideBox(Box<NDIM>::grow(patch->getBox(), fill_ghost_width), axis);
                for (Box<NDIM>::Iterator b(u_data->getArrayData(axis).getBox()); b; b++)
                {
                    if (is_beyond_boundary(axis, b())) continue;
                    double phase = 4.1 * axis;
                    for (unsigned int d = 0; d < NDIM; ++d) phase += (12.9898 + 30.0 * d) * b()(d);
                    u_data->getArrayData(axis)(b(), 0) = fill_box.contains(b()) ?
                                                             std::sin(phase) :
                                                             (run == 0 ? 1.0e3 : -2.5e5) * (1.5 + std::sin(phase));
                }
            }
            bdry_op.setPhysicalBoundaryConditions(*patch, fill_time, IntVector<NDIM>(fill_ghost_width));
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const ArrayData<NDIM, double>& array = u_data->getArrayData(axis);
                const Box<NDIM> fill_box =
                    SideGeometry<NDIM>::toSideBox(Box<NDIM>::grow(patch->getBox(), fill_ghost_width), axis);
                for (Box<NDIM>::Iterator b(array.getBox()); b; b++)
                {
                    // The unfilled faces inside the domain contain values that differ between the runs.
                    if (!is_beyond_boundary(axis, b()) && !fill_box.contains(b())) continue;
                    const double value = array(b(), 0);
                    results[run].push_back(value);
                    if (run == 0 && is_beyond_boundary(axis, b()) && fill_box.contains(b()))
                    {
                        num_finite_outside += std::isfinite(value) ? 1 : 0;
                        num_nan_in_fill_box += std::isnan(value) ? 1 : 0;
                    }
                }
            }
        }
        for (std::size_t n = 0; n < results[0].size(); ++n)
        {
            const bool same =
                (std::isnan(results[0][n]) && std::isnan(results[1][n])) || results[0][n] == results[1][n];
            num_different += same ? 0 : 1;
        }
    }
    level->deallocatePatchData(u_idx);
    if (num_different != 0 || num_finite_outside == 0)
    {
        TBOX_ERROR("check_unfilled_data_unread():\n"
                   << "  " << num_different << " values depend on the faces that the fill does not fill\n"
                   << "  " << num_finite_outside << " values outside the domain are finite\n"
                   << "  the first count must be zero\n"
                   << "  the second count must be positive\n");
    }
    pout << "values that depend on the unfilled faces = " << num_different << '\n'
         << "NaN values outside the domain in the fill boxes = " << num_nan_in_fill_box << '\n';
    return;
} // check_unfilled_data_unread
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // With variable_viscosity, the test uses the viscosity data of the variable-coefficient integrator.
        const bool variable_viscosity = input_db->getBoolWithDefault("variable_viscosity", false);
        Pointer<INSHierarchyIntegrator> ins_integrator;
        if (variable_viscosity)
        {
            ins_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else
        {
            ins_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ins_integrator,
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

        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_coefs[d] =
                new muParserRobinBcCoefs("u_bc_coefs_" + std::to_string(d),
                                         app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                                         grid_geometry);
        }
        ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        Pointer<CartGridFunction> divergence_fcn;
        if (input_db->keyExists("VelocityDivergence"))
        {
            divergence_fcn = new muParserCartGridFunction(
                "divergence_fcn", app_initializer->getComponentDatabase("VelocityDivergence"), grid_geometry);
            ins_integrator->registerVelocityDivergenceFunction(divergence_fcn);
        }

        std::vector<std::unique_ptr<RobinBcCoefStrategy<NDIM>>> mu_bc_coefs;
        if (variable_viscosity)
        {
            Pointer<INSVCStaggeredHierarchyIntegrator> vc_integrator = ins_integrator;
            Pointer<CellVariable<NDIM, double>> mu_var = new CellVariable<NDIM, double>("mu");
            Pointer<CartGridFunction> mu_fcn = new muParserCartGridFunction(
                "mu_fcn", app_initializer->getComponentDatabase("ViscosityFunction"), grid_geometry);
            vc_integrator->registerViscosityVariable(mu_var);
            vc_integrator->registerViscosityInitialConditions(mu_fcn);
            mu_bc_coefs.emplace_back(new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBoundaryConditions"), grid_geometry));
            vc_integrator->registerViscosityBoundaryConditions(mu_bc_coefs.back().get());
        }

        // The filled velocity has ghost width G. Registering it before the integrator, as an IB method does, gives the
        // viscosity data at least that ghost width.
        const int extension_ghost_width = input_db->getInteger("extension_ghost_width");
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
        const int u_idx = var_db->registerVariableAndContext(
            u_var, var_db->getContext("divergence_free_extension"), IntVector<NDIM>(extension_ghost_width));
        const int y_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);

        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // The variable-coefficient integrator sets the viscosity that the operator reads at the start of a time step.
        if (variable_viscosity)
        {
            const double current_time = ins_integrator->getIntegratorTime();
            ins_integrator->preprocessIntegrateHierarchy(
                current_time, current_time + ins_integrator->getMaximumTimeStepSize(), 1);
        }
        check_divergence_free_extension(patch_hierarchy,
                                        ins_integrator,
                                        app_initializer->getComponentDatabase("ExactVelocity"),
                                        extension_ghost_width,
                                        u_idx,
                                        y_idx);
        if (input_db->keyExists("narrow_fill_ghost_width"))
        {
            check_unfilled_data_unread(
                patch_hierarchy, ins_integrator, input_db->getInteger("narrow_fill_ghost_width"), u_idx);
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
