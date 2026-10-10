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

#include <ibamr/INSStaggeredDivergenceFreePhysBdryOp.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
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

    // Evaluate component axis at the position x.
    double operator()(const unsigned int axis, const std::array<double, NDIM>& x)
    {
        d_X = x;
        return d_parsers[axis].Eval();
    }

private:
    std::array<double, NDIM> d_X;
    std::vector<mu::Parser> d_parsers;
};

// Return, for each boundary, whether the velocity condition on it is a traction condition (b = 1) and not a
// prescribed velocity (a = 1). The conditions are constant along each boundary and the same for every component.
std::array<bool, 2 * NDIM>
get_traction_boundaries(Pointer<Database> bc_coefs_db)
{
    std::array<bool, 2 * NDIM> is_traction;
    for (unsigned int location_index = 0; location_index < 2 * NDIM; ++location_index)
    {
        mu::Parser parser;
        parser.SetExpr(bc_coefs_db->getString("bcoef_function_" + std::to_string(location_index)));
        is_traction[location_index] = IBTK::rel_equal_eps(parser.Eval(), 1.0);
    }
    return is_traction;
} // get_traction_boundaries

void
check_divergence_free_extension(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                Pointer<INSStaggeredHierarchyIntegrator> ins_integrator,
                                Pointer<Database> exact_velocity_db,
                                const std::array<bool, 2 * NDIM>& is_traction_boundary,
                                const int extension_ghost_width)
{
    // Fill the ghost values outside the physical boundaries of a level of several patches with the divergence-free
    // extension. Check that the divergence in the ghost cells outside the domain vanishes, that the patches agree on
    // the ghost values, that the values outside the domain are zero after the transpose, and that no ghost value
    // within width G is NaN. Check and report the transpose identity (y, E u) = (E^T y, u), and report the error for a
    // smooth divergence-free velocity whose boundary data match it, separately for the ghost values that lie beyond
    // only boundaries that prescribe the velocity and for those beyond only traction boundaries; the values beyond
    // both kinds are determined by averaging and are checked only through the divergence and the patch agreement. The
    // velocity has ghost width G + 1 so that the ghost values within width G are computable.
    if (IBTK_MPI::getNodes() > 1 || patch_hierarchy->getFinestLevelNumber() > 0)
    {
        TBOX_ERROR(
            "check_divergence_free_extension: the comparison of the values stored by different patches "
            "requires one rank and one level, but there are "
            << IBTK_MPI::getNodes() << " ranks and " << patch_hierarchy->getFinestLevelNumber() + 1 << " levels.\n");
    }
    const int G = extension_ghost_width;
    ParsedVelocity exact_velocity(exact_velocity_db);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const IntVector<NDIM> ghost_width(G + 1);
    const int u_idx =
        var_db->registerVariableAndContext(u_var, var_db->getContext("divergence_free_extension"), ghost_width);
    const int y_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
    INSStaggeredDivergenceFreePhysBdryOp bdry_op(ins_integrator, /*homogeneous_bc*/ false);
    INSStaggeredDivergenceFreePhysBdryOp homogeneous_bdry_op(ins_integrator, /*homogeneous_bc*/ true);
    const double fill_time = ins_integrator->getIntegratorTime();
    const double nan = std::numeric_limits<double>::quiet_NaN();

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
    double prescribed_error = 0.0;
    double traction_error = 0.0;
    double num_prescribed_faces = 0.0;
    double num_traction_faces = 0.0;
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

        // Whether the face i of component axis lies beyond a boundary of the
        // domain. Taking axis = NDIM gives the cell i.
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

        // Call fcn(patch, axis, i) for each side index i of the patch data
        // within width of the patch that is outside the domain.
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

        // Set the values in the patch interiors, and NaN elsewhere so that
        // any dependence on unset values appears in the results.
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

        // Update the maximum divergence in the ghost cells outside the domain, and the scale max |u| / dx of the
        // divergence.
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
                    double divergence = 0.0;
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        hier::Index<NDIM> i_upper = b();
                        i_upper(axis) += 1;
                        const double u_upper = u_data->getArrayData(axis)(i_upper, 0);
                        const double u_lower = u_data->getArrayData(axis)(b(), 0);
                        divergence += (u_upper - u_lower) / pgeom->getDx()[axis];
                        divergence_scale = std::max(
                            divergence_scale, std::max(std::abs(u_upper), std::abs(u_lower)) / pgeom->getDx()[axis]);
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

        // Smooth velocity and boundary data.
        const auto exact = [&](const Patch<NDIM>& patch, const unsigned int axis, const hier::Index<NDIM>& i)
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
        fill_ghosts(exact, ghost_fill_sched);
        update_exterior_divergence();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            for_each_exterior_face(
                patch,
                G,
                [&](const unsigned int axis, const hier::Index<NDIM>& i)
                {
                    const double value = u_data->getArrayData(axis)(i, 0);
                    if (std::isnan(value))
                    {
                        num_nan += 1.0;
                        return;
                    }
                    std::array<int, NDIM + 1> key;
                    key[0] = static_cast<int>(axis);
                    bool beyond_prescribed = false;
                    bool beyond_traction = false;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        key[d + 1] = i(d);
                        if (periodic_shift(d) == 0)
                        {
                            if (i(d) < domain_box.lower(d))
                            {
                                (is_traction_boundary[2 * d] ? beyond_traction : beyond_prescribed) = true;
                            }
                            else if (i(d) > domain_box.upper(d) + (d == axis ? 1 : 0))
                            {
                                (is_traction_boundary[2 * d + 1] ? beyond_traction : beyond_prescribed) = true;
                            }
                        }
                    }
                    const auto it = face_values.find(key);
                    ghost_scale = std::max(ghost_scale, std::abs(value));
                    if (it == face_values.end())
                    {
                        face_values[key] = value;
                    }
                    else
                    {
                        patch_difference = std::max(patch_difference, std::abs(value - it->second));
                    }
                    const double error = std::abs(value - exact(*patch, axis, i));
                    if (beyond_prescribed && !beyond_traction)
                    {
                        num_prescribed_faces += 1.0;
                        prescribed_error = std::max(prescribed_error, error);
                    }
                    else if (beyond_traction && !beyond_prescribed)
                    {
                        num_traction_faces += 1.0;
                        traction_error = std::max(traction_error, error);
                    }
                });
        }

        // Velocity that is not smooth.
        fill_ghosts(rough_value, ghost_fill_sched);
        update_exterior_divergence();

        // Transpose identity. The values y outside the domain within width G are arbitrary, and the values u in
        // the domain are not smooth.
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
                                       y_data->getArrayData(axis)(i, 0) = rough(axis, i, 1.3) + 0.2;
                                       y_dot_Eu += y_data->getArrayData(axis)(i, 0) * u_data->getArrayData(axis)(i, 0);
                                   });
            homogeneous_bdry_op.setPatchDataIndex(y_idx);
            homogeneous_bdry_op.accumulateFromPhysicalBoundaryData(*patch, fill_time, ghost_width);
            homogeneous_bdry_op.setPatchDataIndex(u_idx);

            // The accumulation sets to zero every value outside the domain that
            // the extension fills, which here is every computable value within
            // the ghost width of the patch data.
            for_each_exterior_face(
                patch,
                G + 1,
                [&](const unsigned int axis, const hier::Index<NDIM>& i)
                {
                    const double value = y_data->getArrayData(axis)(i, 0);
                    if (!std::isfinite(value))
                    {
                        TBOX_ERROR("check_divergence_free_extension: the value after the accumulation of component "
                                   << axis << " at the face " << i << " outside the domain is not finite.\n");
                    }
                    max_filled_after_transpose = std::max(max_filled_after_transpose, std::abs(value));
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
        TBOX_ERROR("check_divergence_free_extension: the divergence in a ghost cell outside the domain is "
                   << divergence << ", which exceeds " << divergence_bound << ".\n");
    }
    if (IBTK_MPI::sumReduction(num_nonfinite_divergence) != 0.0 || IBTK_MPI::sumReduction(num_nan) != 0.0)
    {
        TBOX_ERROR("check_divergence_free_extension: a ghost value within width " << G << " is not finite.\n");
    }
    const double patch_difference_max = IBTK_MPI::maxReduction(patch_difference);
    if (!(patch_difference_max <= 1.0e-12 * IBTK_MPI::maxReduction(ghost_scale)))
    {
        TBOX_ERROR("check_divergence_free_extension: two patches differ by " << patch_difference_max
                                                                             << " in a ghost value.\n");
    }
    if (IBTK_MPI::maxReduction(max_filled_after_transpose) != 0.0)
    {
        TBOX_ERROR(
            "check_divergence_free_extension: a value outside the domain is not zero after the "
            "accumulation.\n");
    }
    const double y_dot_Eu_sum = IBTK_MPI::sumReduction(y_dot_Eu);
    const double ETy_dot_u_sum = IBTK_MPI::sumReduction(ETy_dot_u);
    if (!IBTK::rel_equal_eps(y_dot_Eu_sum, ETy_dot_u_sum, 1.0e-12))
    {
        TBOX_ERROR("check_divergence_free_extension: (y, E u) = " << y_dot_Eu_sum << " and (E^T y, u) = "
                                                                  << ETy_dot_u_sum << " do not agree.\n");
    }
    if (IBTK_MPI::sumReduction(num_prescribed_faces) == 0.0 || IBTK_MPI::sumReduction(num_traction_faces) == 0.0)
    {
        TBOX_ERROR(
            "check_divergence_free_extension: the case must have ghost values beyond prescribed-velocity boundaries "
            "and beyond traction boundaries.\n");
    }
    pout << std::setprecision(12) << "(y, E u)   = " << y_dot_Eu_sum << '\n'
         << "(E^T y, u) = " << ETy_dot_u_sum << '\n'
         << "max error beyond prescribed-velocity boundaries = " << IBTK_MPI::maxReduction(prescribed_error) << '\n'
         << "max error beyond traction boundaries            = " << IBTK_MPI::maxReduction(traction_error) << '\n';
    return;
} // check_divergence_free_extension
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<INSStaggeredHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
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

        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // The velocity boundary condition objects of the integrator are set up at the start of a time step.
        const double current_time = ins_integrator->getIntegratorTime();
        ins_integrator->preprocessIntegrateHierarchy(
            current_time, current_time + ins_integrator->getMaximumTimeStepSize(), 1);
        check_divergence_free_extension(
            patch_hierarchy,
            ins_integrator,
            app_initializer->getComponentDatabase("ExactVelocity"),
            get_traction_boundaries(app_initializer->getComponentDatabase("VelocityBcCoefs_0")),
            input_db->getInteger("extension_ghost_width"));

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
