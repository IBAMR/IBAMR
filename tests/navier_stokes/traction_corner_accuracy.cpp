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

#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ArrayData.h>
#include <BergerRigoutsos.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>
#include <muParser.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// The largest error of a set of values and the number of values.
struct RegionError
{
    double max_error = 0.0;
    int count = 0;
};

void
add_error(RegionError& region, const double error)
{
    region.max_error = std::max(region.max_error, error);
    ++region.count;
}

void
check_traction_corner_accuracy(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                               Pointer<INSHierarchyIntegrator> ins_integrator,
                               const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                               const std::string& bdry_interp_type,
                               Pointer<Database> exact_velocity_db)
{
    // Fill the ghost values of a divergence-free velocity field that satisfies the boundary conditions and check the
    // error in the tangential ghost values outside the x boundaries, at corners and elsewhere. At the boundary faces
    // where the normal velocity is not prescribed, check the normal-stress term of the momentum equation and the normal
    // ghost values of the divergence-free extension, next to corners and elsewhere.
    std::array<double, NDIM> X;
    std::vector<mu::Parser> parsers(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        parsers[d].SetExpr(exact_velocity_db->getString("function_" + std::to_string(d)));
        for (unsigned int k = 0; k < NDIM; ++k)
        {
            parsers[d].DefineVar("X_" + std::to_string(k), &X[k]);
        }
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("traction_corner_accuracy");
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const IntVector<NDIM> ghost_width(3);
    const int u_idx = var_db->registerVariableAndContext(u_var, ctx, ghost_width);
    const int u_extension_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
    Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("traction_corner_accuracy_f");
    const int f_idx = var_db->registerVariableAndContext(f_var, ctx, IntVector<NDIM>(0));
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = ins_integrator->getVelocityBoundaryConditions();
    CartSideRobinPhysBdryOp bdry_op(u_idx, U_bc_coefs, /*homogeneous_bc*/ false, bdry_interp_type);
    RobinPhysBdryPatchStrategy* extension_op = ins_integrator->getDivergenceFreeVelocityPhysBdryOp();
    const double fill_time = ins_integrator->getIntegratorTime();
    StaggeredStokesPhysicalBoundaryHelper bc_helper;
    bc_helper.cacheBcCoefData(u_bc_coefs, fill_time, patch_hierarchy);

    RegionError lower_end, upper_end, other, normal_corner, normal_other, stress_corner, stress_other;
    int outermost_count = 0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_idx, fill_time);
        level->allocatePatchData(u_extension_idx, fill_time);
        level->allocatePatchData(f_idx, fill_time);
        BoxArray<NDIM> domain = patch_hierarchy->getGridGeometry()->getPhysicalDomain();
        domain.refine(level->getRatio());
        const Box<NDIM> domain_box = domain.getBoundingBox();
        // Return component axis of the exact velocity at the side i of patch, which may be a ghost side.
        const auto exact_velocity = [&](const Patch<NDIM>& patch, const unsigned int axis, const hier::Index<NDIM>& i)
        {
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X[d] = pgeom->getXLower()[d] +
                       pgeom->getDx()[d] * (i(d) - patch.getBox().lower(d) + (d == axis ? 0.0 : 0.5));
            }
            return parsers[axis].Eval();
        };
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);

            // Set the values in the physical domain, including ghost cells, as a ghost cell fill would.
            u_data->fillAll(0.0);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    hier::Index<NDIM> i_lower = b();
                    i_lower(axis) -= 1;
                    if (domain.contains(b()) || domain.contains(i_lower))
                    {
                        u_data->getArrayData(axis)(b(), 0) = exact_velocity(*patch, axis, b());
                    }
                }
            }
            bdry_op.setPatchDataIndex(u_idx);
            bdry_op.setPhysicalBoundaryConditions(*patch, fill_time, ghost_width);

            // Tangential ghost values outside the x boundaries, including the outermost faces of the patch data.
            for (unsigned int axis = 1; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                const Box<NDIM> ghost_side_box = SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis);
                for (Box<NDIM>::Iterator b(ghost_side_box); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    bool beyond_x_boundary = i(0) < domain_box.lower(0) || i(0) > domain_box.upper(0);
                    bool outermost = false;
                    for (unsigned int d = 1; d < NDIM; ++d)
                    {
                        const bool within_patch = side_box.lower(d) <= i(d) && i(d) <= side_box.upper(d);
                        if (d == axis && !within_patch)
                        {
                            outermost = i(d) == ghost_side_box.lower(d) || i(d) == ghost_side_box.upper(d);
                            beyond_x_boundary = beyond_x_boundary && outermost;
                        }
                        else
                        {
                            beyond_x_boundary = beyond_x_boundary && within_patch;
                        }
                    }
                    if (!beyond_x_boundary)
                    {
                        continue;
                    }
                    // The ghost value is at a corner if the layer i(axis) - 1 or i(axis) of the normal velocity lies
                    // outside the physical domain.
                    hier::Index<NDIM> i_boundary = i;
                    i_boundary(0) = i(0) < patch_box.lower(0) ? patch_box.lower(0) : patch_box.upper(0);
                    hier::Index<NDIM> i_boundary_lower = i_boundary;
                    i_boundary_lower(axis) -= 1;
                    const double error = std::abs(u_data->getArrayData(axis)(i, 0) - exact_velocity(*patch, axis, i));
                    const bool lower_layer_in_domain = domain.contains(i_boundary_lower);
                    const bool upper_layer_in_domain = domain.contains(i_boundary);
                    if (outermost)
                    {
                        // The fill sets the outermost ghost faces to NaN, because their difference involves a normal
                        // velocity face outside the patch data.
                        if (lower_layer_in_domain && upper_layer_in_domain)
                        {
                            if (!std::isnan(u_data->getArrayData(axis)(i, 0)))
                            {
                                TBOX_ERROR("check_traction_corner_accuracy():\n"
                                           << "  the outermost ghost face " << i << " of axis " << axis
                                           << " is not NaN\n");
                            }
                            ++outermost_count;
                        }
                    }
                    else if (!lower_layer_in_domain)
                    {
                        add_error(lower_end, error);
                    }
                    else if (!upper_layer_in_domain)
                    {
                        add_error(upper_end, error);
                    }
                    else
                    {
                        add_error(other, error);
                    }
                }
            }

            // The normal ghost values of the divergence-free extension, which replaces all ghost values outside the
            // physical boundaries.
            Pointer<SideData<NDIM, double>> u_extension_data = patch->getPatchData(u_extension_idx);
            u_extension_data->copy(*u_data);
            extension_op->setPatchDataIndex(u_extension_idx);
            extension_op->setPhysicalBoundaryConditions(*patch, fill_time, IntVector<NDIM>(1));
            Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
            f_data->fillAll(0.0);
        }

        // Add the normal-stress term of the momentum equation to zero data.
        bc_helper.addNormalTractionViscousTerm(
            f_idx, u_idx, /*viscous_coef*/ 1.0, U_bc_coefs, bdry_interp_type, /*coarsest_ln*/ ln, /*finest_ln*/ ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            Pointer<SideData<NDIM, double>> u_extension_data = patch->getPatchData(u_extension_idx);
            Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
            const tbox::Array<BoundaryBox<NDIM>>& bdry_boxes = pgeom->getCodimensionBoundaries(1);
            for (int k = 0; k < bdry_boxes.size(); ++k)
            {
                const unsigned int axis = bdry_boxes[k].getLocationIndex() / 2;
                const bool is_lower = bdry_boxes[k].getLocationIndex() % 2 == 0;
                const double h = pgeom->getDx()[axis];
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_boxes[k]);
                Pointer<ArrayData<NDIM, double>> a_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> b_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> g_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                u_bc_coefs[axis]->setBcCoefs(
                    a_data, b_data, g_data, Pointer<Variable<NDIM>>(), *patch, bdry_boxes[k], fill_time);

                // The cell it() is the interior cell abutting a lower boundary or the ghost cell abutting an upper
                // boundary, so that its lower face is the boundary face.
                Box<NDIM> cell_box = patch->getBox();
                cell_box.lower(axis) -= 1;
                cell_box.upper(axis) += 1;
                for (Box<NDIM>::Iterator it(bc_coef_box * cell_box); it; it++)
                {
                    if (!((*a_data)(it(), 0) == 0.0 && (*b_data)(it(), 0) == 1.0))
                    {
                        continue;
                    }
                    // The face is next to a corner if the interior cell abutting it has a neighbor outside the
                    // physical domain along another axis.
                    hier::Index<NDIM> i_cell = it();
                    if (!is_lower)
                    {
                        i_cell(axis) -= 1;
                    }
                    bool at_corner = false;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        hier::Index<NDIM> i_below = i_cell, i_above = i_cell;
                        i_below(d) -= 1;
                        i_above(d) += 1;
                        at_corner =
                            at_corner || (d != axis && (!domain.contains(i_below) || !domain.contains(i_above)));
                    }
                    const int inward = is_lower ? 1 : -1;
                    SideIndex<NDIM> s_boundary(it(), axis, SideIndex<NDIM>::Lower), s_inside(s_boundary),
                        s_inside_2(s_boundary), s_outside(s_boundary);
                    s_inside(axis) += inward;
                    s_inside_2(axis) += 2 * inward;
                    s_outside(axis) -= inward;
                    const double u_boundary = exact_velocity(*patch, axis, s_boundary);
                    const double u_inside = exact_velocity(*patch, axis, s_inside);
                    const double u_inside_2 = exact_velocity(*patch, axis, s_inside_2);
                    const double u_outside = exact_velocity(*patch, axis, s_outside);
                    // The term with the exact values, extrapolated as the pressure ghost values are.
                    const double exact_term = bdry_interp_type == "QUADRATIC" ?
                                                  (3.0 * (u_inside - u_boundary) -
                                                   (2.0 / 3.0) * (u_inside_2 - u_inside) + (u_boundary - u_outside)) /
                                                      (h * h) :
                                                  (u_inside - u_outside) / (h * h);
                    const double u_extension_outside = (*u_extension_data)(s_outside);
                    const double term = (*f_data)(s_boundary);
                    const double normal_error = std::abs(u_extension_outside - u_outside);
                    const double stress_error = std::abs(term - exact_term);
                    add_error(at_corner ? normal_corner : normal_other, normal_error);
                    add_error(at_corner ? stress_corner : stress_other, stress_error);
                }
            }
        }
        level->deallocatePatchData(u_idx);
        level->deallocatePatchData(u_extension_idx);
        level->deallocatePatchData(f_idx);
    }

    // Check the largest errors against the bounds in the input and print the number of values checked and the largest
    // error in each region.
    for (RegionError* region :
         { &lower_end, &upper_end, &other, &normal_corner, &normal_other, &stress_corner, &stress_other })
    {
        region->max_error = IBTK_MPI::maxReduction(region->max_error);
        region->count = IBTK_MPI::sumReduction(region->count);
    }
    const double error_bound = exact_velocity_db->getDouble("error_bound");
    const double max_error = std::max(
        { lower_end.max_error, upper_end.max_error, other.max_error, normal_corner.max_error, normal_other.max_error });
    if (!(max_error <= error_bound))
    {
        TBOX_ERROR("check_traction_corner_accuracy():\n"
                   << "  the largest ghost error " << max_error << " exceeds the bound " << error_bound << "\n");
    }
    const double stress_error_bound = exact_velocity_db->getDouble("stress_error_bound");
    const double max_stress_error = std::max(stress_corner.max_error, stress_other.max_error);
    if (!(max_stress_error <= stress_error_bound))
    {
        TBOX_ERROR("check_traction_corner_accuracy():\n"
                   << "  the largest normal-stress term error " << max_stress_error << " exceeds the bound "
                   << stress_error_bound << "\n");
    }
    const auto print = [](const std::string& description, const RegionError& region)
    {
        plog << description << ": " << region.count << " checked, largest error " << std::fixed << std::setprecision(10)
             << region.max_error << '\n';
    };
    print("ghost faces at the lower corners of the x boundaries", lower_end);
    print("ghost faces at the upper corners of the x boundaries", upper_end);
    print("ghost faces elsewhere on the x boundaries", other);
    plog << "outermost ghost faces at the ends of the ghost box: " << IBTK_MPI::sumReduction(outermost_count)
         << " are NaN, as the fill does not compute them\n";
    print("normal velocity ghost faces of the extension next to corners", normal_corner);
    print("normal velocity ghost faces of the extension elsewhere", normal_other);
    print("normal-stress terms next to corners", stress_corner);
    print("normal-stress terms elsewhere", stress_other);
    return;
} // check_traction_corner_accuracy
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");

        Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
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

        ins_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // The integrator sets up its velocity boundary condition objects at the start of a time step.
        const double current_time = ins_integrator->getIntegratorTime();
        ins_integrator->preprocessIntegrateHierarchy(
            current_time, current_time + ins_integrator->getMaximumTimeStepSize(), 1);
        const std::string bdry_interp_type = app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator")
                                                 ->getStringWithDefault("bdry_interp_type", "LINEAR");
        check_traction_corner_accuracy(patch_hierarchy,
                                       ins_integrator,
                                       u_bc_coefs,
                                       bdry_interp_type,
                                       app_initializer->getComponentDatabase("ExactVelocity"));
    }
} // main
