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
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/StaggeredStokesBoxRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// Smooth, nonzero values of the error.
double
exact_value(const hier::Index<NDIM>& i, const int component)
{
    double arg = 0.8 * component;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        arg += (1.7 + 1.2 * d) * i(d);
    }
    return std::sin(arg);
}

// Return the largest magnitude of the momentum and continuity components of the vector with patch data indices
// (u_idx, p_idx) over the patch interiors.
double
largest_magnitude(Pointer<PatchLevel<NDIM>> level, const int u_idx, const int p_idx)
{
    double result = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
            {
                const double value = (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower));
                if (!std::isfinite(value))
                {
                    TBOX_ERROR("a velocity value is not finite\n");
                }
                result = std::max(result, std::abs(value));
            }
        }
        for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
        {
            const double value = (*p_data)(it());
            if (!std::isfinite(value))
            {
                TBOX_ERROR("a pressure value is not finite\n");
            }
            result = std::max(result, std::abs(value));
        }
    }
    return IBTK_MPI::maxReduction(result);
}

// Take an error e that satisfies the homogeneous boundary conditions, and the residual r = A*e of the staggered Stokes
// operator A with the integrator's boundary condition objects. The residual that the FAC strategy computes for the
// solution e and the right-hand side zero must be -r in every row, and one sweep of the box relaxation smoother for
// the right-hand side r must leave e unchanged.
void
check_fac_open_boundary(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                        Pointer<INSHierarchyIntegrator> ins_integrator,
                        const vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    const double solution_time = 0.5;
    const double mu = ins_integrator->getStokesSpecifications()->getMu();

    // Configure the integrator's boundary condition objects.
    const vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = ins_integrator->getVelocityBoundaryConditions();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto U_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(U_bc_coefs[d]);
        U_bc_coef->setStokesSpecifications(ins_integrator->getStokesSpecifications());
        U_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
        U_bc_coef->setSolutionTime(solution_time);
    }
    RobinBcCoefStrategy<NDIM>* P_bc_coef = ins_integrator->getPressureBoundaryConditions();
    auto P_ins_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(P_bc_coef);
    P_ins_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
    P_ins_bc_coef->setSolutionTime(solution_time);

    // The error e, its residual r, the residual that the FAC strategy computes, the zero right-hand side, and a copy of
    // the error.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("fac_open_boundary");
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<SideVariable<NDIM, double>>> u_vars;
    std::vector<Pointer<CellVariable<NDIM, double>>> p_vars;
    std::vector<int> u_idxs, p_idxs;
    for (const char* name : { "e", "r", "r_fac", "zero", "e_copy", "difference" })
    {
        u_vars.push_back(new SideVariable<NDIM, double>(std::string("u_") + name));
        p_vars.push_back(new CellVariable<NDIM, double>(std::string("p_") + name));
        u_idxs.push_back(var_db->registerVariableAndContext(u_vars.back(), ctx, IntVector<NDIM>(1)));
        p_idxs.push_back(var_db->registerVariableAndContext(p_vars.back(), ctx, IntVector<NDIM>(1)));
    }
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> vecs;
    for (std::size_t k = 0; k < u_vars.size(); ++k)
    {
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_idxs[k], solution_time);
            level->allocatePatchData(p_idxs[k], solution_time);
        }
        vecs.push_back(new SAMRAIVectorReal<NDIM, double>("vec_" + std::to_string(k), patch_hierarchy, 0, finest_ln));
        vecs.back()->addComponent(u_vars[k], u_idxs[k]);
        vecs.back()->addComponent(p_vars[k], p_idxs[k]);
        vecs.back()->setToScalar(0.0, /*interior_only*/ false);
    }
    SAMRAIVectorReal<NDIM, double>& e = *vecs[0];
    SAMRAIVectorReal<NDIM, double>& r = *vecs[1];
    SAMRAIVectorReal<NDIM, double>& r_fac = *vecs[2];
    SAMRAIVectorReal<NDIM, double>& zero = *vecs[3];
    SAMRAIVectorReal<NDIM, double>& e_copy = *vecs[4];
    SAMRAIVectorReal<NDIM, double>& difference = *vecs[5];
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idxs[0]);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idxs[0]);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
                {
                    (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) = exact_value(it(), axis);
                }
            }
            for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
            {
                (*p_data)(it()) = exact_value(it(), NDIM);
            }
        }
    }

    // The error satisfies the homogeneous boundary conditions: the normal velocity vanishes where it is prescribed.
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    bc_helper->cacheBcCoefData(u_bc_coefs, solution_time, patch_hierarchy);
    bc_helper->enforceNormalVelocityBoundaryConditions(
        u_idxs[0], p_idxs[0], U_bc_coefs, solution_time, /*homogeneous_bc*/ true);

    // The residual is the operator applied to the error with homogeneous boundary conditions.
    PoissonSpecifications spec("spec");
    spec.setCConstant(1.0);
    spec.setDConstant(-mu);
    StaggeredStokesOperator op("StaggeredStokesOperator", /*homogeneous_bc*/ true);
    op.setVelocityPoissonSpecifications(spec);
    op.setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    op.setPhysicalBoundaryHelper(bc_helper);
    op.setSolutionTime(solution_time);
    op.setTimeInterval(solution_time, solution_time);
    op.initializeOperatorState(e, r);
    op.apply(e, r);
    op.deallocateOperatorState();

    StaggeredStokesBoxRelaxationFACOperator fac_op("box_smoother", nullptr, "box_");
    fac_op.setVelocityPoissonSpecifications(spec);
    fac_op.setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    fac_op.setPhysicalBoundaryHelper(bc_helper);
    fac_op.setSolutionTime(solution_time);
    fac_op.setTimeInterval(solution_time, solution_time);
    fac_op.initializeOperatorState(e, r);

    // The FAC strategy's residual for the right-hand side zero is minus the operator's residual.
    fac_op.computeResidual(r_fac, e, zero, 0, finest_ln);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    difference.add(Pointer<SAMRAIVectorReal<NDIM, double>>(&r_fac, false),
                   Pointer<SAMRAIVectorReal<NDIM, double>>(&r, false));
    const double r_norm = largest_magnitude(level, u_idxs[1], p_idxs[1]);
    const double residual_difference = largest_magnitude(level, u_idxs[5], p_idxs[5]);
    if (residual_difference > 1.0e-10 * r_norm)
    {
        TBOX_ERROR("the residual of the FAC strategy differs from the residual of the operator by "
                   << residual_difference << "; the largest magnitude of the operator's residual is " << r_norm
                   << "\n");
    }

    // The error is a fixed point of the smoother for the right-hand side r. The sweep changes the error in the patch
    // interiors only, and its ghost values were filled by the residual computation.
    const double e_norm = largest_magnitude(level, u_idxs[0], p_idxs[0]);
    e_copy.copyVector(Pointer<SAMRAIVectorReal<NDIM, double>>(&e, false), /*interior_only*/ false);
    fac_op.smoothError(e, r, 0, 1, false, false);
    difference.subtract(Pointer<SAMRAIVectorReal<NDIM, double>>(&e, false),
                        Pointer<SAMRAIVectorReal<NDIM, double>>(&e_copy, false));
    const double change = largest_magnitude(level, u_idxs[5], p_idxs[5]);
    if (change > 1.0e-10 * e_norm)
    {
        TBOX_ERROR("one sweep of the smoother changes the exact error by "
                   << change << "; the largest magnitude of the error is " << e_norm << "\n");
    }
    fac_op.deallocateOperatorState();
    plog << std::setprecision(10) << "Largest magnitude of the residual: " << r_norm << "\n";
} // check_fac_open_boundary
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
        check_fac_open_boundary(patch_hierarchy, ins_integrator, u_bc_coefs);

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
