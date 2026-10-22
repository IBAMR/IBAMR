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

#include <tbox/DatabaseBox.h>

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
#include <memory>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// The velocity given by the muParser expressions function_0, function_1, ... of the position X_0, X_1, ....
class ExactVelocity
{
public:
    explicit ExactVelocity(Pointer<Database> db) : d_parsers(NDIM)
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

    ExactVelocity(const ExactVelocity&) = delete;
    ExactVelocity& operator=(const ExactVelocity&) = delete;

    // Return component axis at the side i of patch, which may be a ghost side.
    double operator()(const Patch<NDIM>& patch, const unsigned int axis, const hier::Index<NDIM>& i)
    {
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_X[d] =
                pgeom->getXLower()[d] + pgeom->getDx()[d] * (i(d) - patch.getBox().lower(d) + (d == axis ? 0.0 : 0.5));
        }
        return d_parsers[axis].Eval();
    }

private:
    std::array<double, NDIM> d_X;
    std::vector<mu::Parser> d_parsers;
};

// The largest errors in the ghost values and the number of values checked. The tangential velocity ghost values
// outside the x boundaries are checked at the lower and upper corners of these boundaries, on the outermost ghost
// faces of the patch data, and elsewhere. The normal velocity ghost values and the normal-stress term of the momentum
// equation are checked at the rows of all boundaries next to a corner and elsewhere.
struct Errors
{
    double lower_end = 0.0;
    double upper_end = 0.0;
    double other = 0.0;
    double outermost = 0.0;
    double normal_corner = 0.0;
    double normal_other = 0.0;
    double stress_corner = 0.0;
    double stress_other = 0.0;
    int lower_end_count = 0;
    int upper_end_count = 0;
    int other_count = 0;
    int outermost_count = 0;
    int normal_corner_count = 0;
    int normal_other_count = 0;
};

Errors
check_traction_corner_accuracy(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                               Pointer<INSHierarchyIntegrator> ins_integrator,
                               const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                               ExactVelocity& exact_velocity)
{
    // Fill the ghost values of a velocity field that satisfies the boundary conditions and find the errors in the
    // tangential ghost values outside the x boundaries, at corners and elsewhere, and in the normal velocity ghost
    // values and the normal-stress term next to the boundaries.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const IntVector<NDIM> ghost_width(3);
    const int u_idx =
        var_db->registerVariableAndContext(u_var, var_db->getContext("traction_corner_accuracy"), ghost_width);
    Pointer<Variable<NDIM>> f_var = var_db->getVariable("traction_corner_accuracy_f");
    if (!f_var)
    {
        f_var = new SideVariable<NDIM, double>("traction_corner_accuracy_f");
    }
    const int f_idx =
        var_db->registerVariableAndContext(f_var, var_db->getContext("traction_corner_accuracy"), IntVector<NDIM>(0));
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = ins_integrator->getVelocityBoundaryConditions();
    CartSideRobinPhysBdryOp bdry_op(u_idx, U_bc_coefs, /*homogeneous_bc*/ false);
    const double fill_time = ins_integrator->getIntegratorTime();
    StaggeredStokesPhysicalBoundaryHelper bc_helper;
    bc_helper.cacheBcCoefData(u_bc_coefs, fill_time, patch_hierarchy);

    Errors errors;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_idx, fill_time);
        level->allocatePatchData(f_idx, fill_time);
        BoxArray<NDIM> domain = patch_hierarchy->getGridGeometry()->getPhysicalDomain();
        domain.refine(level->getRatio());
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            const auto exact = [&](const unsigned int axis, const hier::Index<NDIM>& i)
            { return exact_velocity(*patch, axis, i); };

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
                        u_data->getArrayData(axis)(b(), 0) = exact(axis, b());
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
                    bool beyond_x_boundary = i(0) < patch_box.lower(0) || i(0) > patch_box.upper(0);
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
                    // The ghost value uses the normal velocity in the rows i(axis) - 1 and i(axis); it is at a
                    // corner if one of them is outside the physical domain.
                    hier::Index<NDIM> i_boundary = i;
                    i_boundary(0) = i(0) < patch_box.lower(0) ? patch_box.lower(0) : patch_box.upper(0);
                    hier::Index<NDIM> i_boundary_lower = i_boundary;
                    i_boundary_lower(axis) -= 1;
                    const double error = std::abs(u_data->getArrayData(axis)(i, 0) - exact(axis, i));
                    const bool lower_row_in_domain = domain.contains(i_boundary_lower);
                    const bool upper_row_in_domain = domain.contains(i_boundary);
                    if (outermost)
                    {
                        if (lower_row_in_domain && upper_row_in_domain)
                        {
                            errors.outermost = std::max(errors.outermost, error);
                            ++errors.outermost_count;
                        }
                    }
                    else if (!lower_row_in_domain)
                    {
                        errors.lower_end = std::max(errors.lower_end, error);
                        ++errors.lower_end_count;
                    }
                    else if (!upper_row_in_domain)
                    {
                        errors.upper_end = std::max(errors.upper_end, error);
                        ++errors.upper_end_count;
                    }
                    else
                    {
                        errors.other = std::max(errors.other, error);
                        ++errors.other_count;
                    }
                }
            }
        }

        // Set the normal velocity ghost values next to the TRACTION boundaries so that the divergence vanishes in
        // the ghost cells, and compute the normal-stress term on a zero field.
        bc_helper.enforceDivergenceFreeConditionAtBoundary(u_idx, ln, ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<SideData<NDIM, double>> f_data = level->getPatch(p())->getPatchData(f_idx);
            f_data->fillAll(0.0);
        }
        bc_helper.addNormalTractionViscousTerm(
            f_idx, u_idx, /*viscous_coef*/ 1.0, U_bc_coefs, /*linear_pressure_extrapolation*/ true, ln, ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
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

                // The boundary faces of the patch. The face of a lower boundary is the lower face of the interior
                // cell abutting it; the face of an upper boundary is the lower face of the ghost cell abutting it.
                Box<NDIM> cell_box = patch->getBox();
                cell_box.lower(axis) -= 1;
                cell_box.upper(axis) += 1;
                for (Box<NDIM>::Iterator it(bc_coef_box * cell_box); it; it++)
                {
                    if (!((*a_data)(it(), 0) == 0.0 && (*b_data)(it(), 0) == 1.0))
                    {
                        continue;
                    }
                    // A row is next to a corner if the cell abutting the boundary at its face has a neighbor along a
                    // boundary that is outside the physical domain.
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
                    SideIndex<NDIM> s_boundary(it(), axis, SideIndex<NDIM>::Lower), s_inside(s_boundary),
                        s_outside(s_boundary);
                    s_inside(axis) += is_lower ? 1 : -1;
                    s_outside(axis) += is_lower ? -1 : 1;
                    const double normal_error = std::abs((*u_data)(s_outside)-exact_velocity(*patch, axis, s_outside));
                    const double stress_error = std::abs(
                        (*f_data)(s_boundary) -
                        (exact_velocity(*patch, axis, s_inside) - exact_velocity(*patch, axis, s_outside)) / (h * h));
                    if (at_corner)
                    {
                        errors.normal_corner = std::max(errors.normal_corner, normal_error);
                        errors.stress_corner = std::max(errors.stress_corner, stress_error);
                        ++errors.normal_corner_count;
                    }
                    else
                    {
                        errors.normal_other = std::max(errors.normal_other, normal_error);
                        errors.stress_other = std::max(errors.stress_other, stress_error);
                        ++errors.normal_other_count;
                    }
                }
            }
        }
        level->deallocatePatchData(u_idx);
        level->deallocatePatchData(f_idx);
    }

    // Combine the values of all ranks.
    errors.lower_end = IBTK_MPI::maxReduction(errors.lower_end);
    errors.upper_end = IBTK_MPI::maxReduction(errors.upper_end);
    errors.other = IBTK_MPI::maxReduction(errors.other);
    errors.outermost = IBTK_MPI::maxReduction(errors.outermost);
    errors.normal_corner = IBTK_MPI::maxReduction(errors.normal_corner);
    errors.normal_other = IBTK_MPI::maxReduction(errors.normal_other);
    errors.stress_corner = IBTK_MPI::maxReduction(errors.stress_corner);
    errors.stress_other = IBTK_MPI::maxReduction(errors.stress_other);
    errors.lower_end_count = IBTK_MPI::sumReduction(errors.lower_end_count);
    errors.upper_end_count = IBTK_MPI::sumReduction(errors.upper_end_count);
    errors.other_count = IBTK_MPI::sumReduction(errors.other_count);
    errors.outermost_count = IBTK_MPI::sumReduction(errors.outermost_count);
    errors.normal_corner_count = IBTK_MPI::sumReduction(errors.normal_corner_count);
    errors.normal_other_count = IBTK_MPI::sumReduction(errors.normal_other_count);
    return errors;
} // check_traction_corner_accuracy

// Build the hierarchy of the input and find the errors on it.
Errors
run(Pointer<AppInitializer> app_initializer, ExactVelocity& exact_velocity)
{
    Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", ins_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
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
    return check_traction_corner_accuracy(patch_hierarchy, ins_integrator, u_bc_coefs, exact_velocity);
} // run
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> exact_velocity_db = app_initializer->getComponentDatabase("ExactVelocity");
        ExactVelocity exact_velocity(exact_velocity_db);
        const Errors errors = run(app_initializer, exact_velocity);

        // Check the largest errors against the bounds in the input and print the number of faces checked.
        const double error_bound = exact_velocity_db->getDouble("error_bound");
        const double max_error =
            std::max(std::max(std::max(errors.lower_end, errors.upper_end), std::max(errors.other, errors.outermost)),
                     std::max(errors.normal_corner, errors.normal_other));
        if (!(max_error <= error_bound))
        {
            TBOX_ERROR("traction_corner_accuracy: the largest ghost error " << max_error << " exceeds the bound "
                                                                            << error_bound << "\n");
        }
        const double stress_error_bound = exact_velocity_db->getDouble("stress_error_bound");
        const double max_stress_error = std::max(errors.stress_corner, errors.stress_other);
        if (!(max_stress_error <= stress_error_bound))
        {
            TBOX_ERROR("traction_corner_accuracy: the largest normal-stress term error "
                       << max_stress_error << " exceeds the bound " << stress_error_bound << "\n");
        }
        pout << "ghost faces checked at the lower corners of the x boundaries = " << errors.lower_end_count << '\n'
             << "ghost faces checked at the upper corners of the x boundaries = " << errors.upper_end_count << '\n'
             << "ghost faces checked elsewhere on the x boundaries            = " << errors.other_count << '\n'
             << "ghost faces checked at the outermost ghost faces             = " << errors.outermost_count << '\n'
             << "normal velocity ghost faces checked at the corner rows       = " << errors.normal_corner_count << '\n'
             << "normal velocity ghost faces checked elsewhere                = " << errors.normal_other_count << '\n';
    }
} // main
