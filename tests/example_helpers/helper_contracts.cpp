// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/PIO.h>

#include <BergerRigoutsos.h>
#include <CapillaryForces.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <GriddingAlgorithm.h>
#include <HeavisideFromLevelSet.h>
#include <HierarchyCellDataOpsReal.h>
#include <LiquidFractionForceMask.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
void
check_value(const double actual, const double expected)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1.0e-12)
    {
        TBOX_ERROR("Helper value " << actual << " differs from expected " << expected << '\n');
    }
}

void
check_heaviside(Pointer<PatchHierarchy<NDIM>> hierarchy,
                Pointer<AdvDiffHierarchyIntegrator> integrator,
                Pointer<CellVariable<NDIM, double>> ls_var,
                const int H_idx)
{
    VariableDatabase<NDIM>* db = VariableDatabase<NDIM>::getDatabase();
    const int ls_idx = db->mapVariableAndContextToIndex(ls_var, integrator->getCurrentContext());
    const std::array<double, 5> distances = { -2.0, -0.5, 0.0, 0.5, 2.0 };
    const double pi = std::acos(-1.0);
    const std::array<double, 5> expected = { 0.0, 0.25 - 0.5 / pi, 0.5, 0.75 + 0.5 / pi, 1.0 };
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
        double volume = 1.0;
        for (int d = 0; d < NDIM; ++d)
        {
            volume *= geometry->getDx()[d];
        }
        const double alpha = 2.0 * std::pow(volume, 1.0 / NDIM);
        Pointer<CellData<NDIM, double>> ls_data = patch->getPatchData(ls_idx);
        for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
        {
            (*ls_data)(CellIndex<NDIM>(i())) = alpha * distances[i()(0) % distances.size()];
        }
    }
    PhaseChangeExamples::HeavisideFromLevelSet context(integrator, ls_var, 2.0);
    Pointer<HierarchyMathOps> math_ops = new HierarchyMathOps("heaviside_math", hierarchy);
    PhaseChangeExamples::HeavisideFromLevelSet::synchronize_levelset_with_heaviside_fcn(
        H_idx, math_ops, 0, 0.0, true, false, &context);
    std::array<double, 5> measured = {};
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_idx);
        for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
        {
            const auto sample = i()(0) % expected.size();
            measured[sample] = (*H_data)(CellIndex<NDIM>(i()));
            check_value(measured[sample], expected[sample]);
        }
    }
    plog << "Heaviside:";
    for (const double value : measured)
    {
        plog << ' ' << value;
    }
    plog << '\n';
}

enum class CapillaryOperation
{
    SURFACE_TENSION,
    MARANGONI,
    DENSITY_MASK,
    LIQUID_MASK,
    EXTRAPOLATED_MASK
};

void
check_capillary(Pointer<PatchHierarchy<NDIM>> hierarchy,
                Pointer<INSVCStaggeredHierarchyIntegrator> ins_integrator,
                Pointer<AdvDiffHierarchyIntegrator> integrator,
                Pointer<CellVariable<NDIM, double>> T_var,
                Pointer<CellVariable<NDIM, double>> lf_var,
                Pointer<CellVariable<NDIM, double>> extrap_var)
{
    VariableDatabase<NDIM>* db = VariableDatabase<NDIM>::getDatabase();
    const int T_idx = db->mapVariableAndContextToIndex(T_var, integrator->getScratchContext());
    const int lf_new_idx = db->mapVariableAndContextToIndex(lf_var, integrator->getNewContext());
    const int lf_scratch_idx = db->mapVariableAndContextToIndex(lf_var, integrator->getScratchContext());
    const int extrap_new_idx = db->mapVariableAndContextToIndex(extrap_var, integrator->getNewContext());
    const int extrap_scratch_idx = db->mapVariableAndContextToIndex(extrap_var, integrator->getScratchContext());
    const int rho_idx = ins_integrator->getLinearOperatorRhoPatchDataIndex();
    Pointer<SideVariable<NDIM, double>> force_var = new SideVariable<NDIM, double>("contract_force");
    const int force_idx = db->registerVariableAndContext(force_var, db->getContext("contract"));
    const std::vector<int> temporary_indices = { T_idx,          lf_new_idx,         lf_scratch_idx,
                                                 extrap_new_idx, extrap_scratch_idx, rho_idx,
                                                 force_idx };
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    for (const int idx : temporary_indices)
    {
        TBOX_ASSERT(!level->checkAllocated(idx));
        level->allocatePatchData(idx, 0.0);
    }
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> T_data = patch->getPatchData(T_idx);
        for (Box<NDIM>::Iterator i(T_data->getGhostBox()); i; i++)
        {
            (*T_data)(CellIndex<NDIM>(i())) = 1.0 + 0.25 * i()(0);
        }
        Pointer<CellData<NDIM, double>> lf_data = patch->getPatchData(lf_new_idx);
        Pointer<CellData<NDIM, double>> extrap_data = patch->getPatchData(extrap_new_idx);
        lf_data->fillAll(0.25);
        extrap_data->fillAll(0.75);
        Pointer<SideData<NDIM, double>> rho_data = patch->getPatchData(rho_idx);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator i(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); i; i++)
            {
                (*rho_data)(SideIndex<NDIM>(i(), axis, SideIndex<NDIM>::Lower)) = 2.0 + axis;
            }
        }
    }
    MultiphaseExamples::SurfaceTensionCoefficients coefficients(T_var, integrator->getScratchContext(), 2.0, -0.5, 1.0);
    MultiphaseExamples::DensityForceMask density_mask(ins_integrator, 3.0, 1.0);
    PhaseChangeExamples::LiquidFractionForceMask liquid_mask(lf_var, nullptr, integrator, ins_integrator, 3.0, 1.0);
    PhaseChangeExamples::LiquidFractionForceMask extrapolated_mask(
        extrap_var, nullptr, integrator, ins_integrator, 3.0, 1.0);
    Pointer<HierarchyMathOps> math_ops = new HierarchyMathOps("capillary_math", hierarchy);
    for (const auto operation : { CapillaryOperation::SURFACE_TENSION,
                                  CapillaryOperation::MARANGONI,
                                  CapillaryOperation::DENSITY_MASK,
                                  CapillaryOperation::LIQUID_MASK,
                                  CapillaryOperation::EXTRAPOLATED_MASK })
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> force_data = patch->getPatchData(force_idx);
            force_data->fillAll(2.0);
            if (operation == CapillaryOperation::SURFACE_TENSION)
            {
                MultiphaseExamples::SurfaceTensionCoefficients::compute_surface_tension_coef_function(
                    force_idx, patch, 0, 0.0, 0.0, 1.0, &coefficients);
            }
            else if (operation == CapillaryOperation::MARANGONI)
            {
                MultiphaseExamples::SurfaceTensionCoefficients::compute_marangoni_coef_function(
                    force_idx, patch, 0, 0.0, 0.0, 1.0, &coefficients);
            }
        }
        if (operation == CapillaryOperation::DENSITY_MASK)
        {
            MultiphaseExamples::DensityForceMask::mask_surface_tension_force(
                force_idx, math_ops, 0, 0.0, 0.0, 1.0, &density_mask);
        }
        else if (operation == CapillaryOperation::LIQUID_MASK || operation == CapillaryOperation::EXTRAPOLATED_MASK)
        {
            auto* context = operation == CapillaryOperation::LIQUID_MASK ? &liquid_mask : &extrapolated_mask;
            PhaseChangeExamples::LiquidFractionForceMask::mask_surface_tension_force(
                force_idx, math_ops, 0, 0.0, 0.0, 1.0, context);
        }
        std::array<double, NDIM> measured = {};
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> force_data = patch->getPatchData(force_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
                for (Box<NDIM>::Iterator i(side_box); i; i++)
                {
                    const SideIndex<NDIM> side(i(), axis, SideIndex<NDIM>::Lower);
                    double expected = 0.0;
                    switch (operation)
                    {
                    case CapillaryOperation::SURFACE_TENSION:
                        expected = 4.0 - 0.25 * (i()(0) - (axis == 0 ? 0.5 : 0.0));
                        break;
                    case CapillaryOperation::MARANGONI:
                        expected = -1.0;
                        break;
                    case CapillaryOperation::DENSITY_MASK:
                        expected = 2.0 + axis;
                        break;
                    case CapillaryOperation::LIQUID_MASK:
                        expected = (2.0 + axis) * 0.25;
                        break;
                    case CapillaryOperation::EXTRAPOLATED_MASK:
                        expected = (2.0 + axis) * 0.75;
                        break;
                    default:
                        TBOX_ERROR("Unknown capillary operation\n");
                    }
                    check_value((*force_data)(side), expected);
                }
                measured[axis] = (*force_data)(SideIndex<NDIM>(side_box.lower(), axis, SideIndex<NDIM>::Lower));
            }
        }
        const std::array<const char*, 5> names = {
            "Surface tension", "Marangoni", "Density mask", "Liquid mask", "Extrapolated mask"
        };
        plog << names[static_cast<std::size_t>(operation)] << ':';
        for (const double value : measured)
        {
            plog << ' ' << value;
        }
        plog << '\n';
    }
    for (const int idx : temporary_indices)
    {
        level->deallocatePatchData(idx);
    }
    db->removePatchDataIndex(force_idx);
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    {
        TimerManager::createManager(nullptr);
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "helper_setup.log");
        Pointer<Database> input = app->getInputDatabase();
        Pointer<INSVCStaggeredHierarchyIntegrator> ins =
            new INSVCStaggeredConservativeHierarchyIntegrator("ins", app->getComponentDatabase("INS"));
        Pointer<AdvDiffHierarchyIntegrator> adv =
            new AdvDiffSemiImplicitHierarchyIntegrator("adv", app->getComponentDatabase("AdvDiff"));
        ins->registerAdvDiffHierarchyIntegrator(adv);
        Pointer<CartesianGridGeometry<NDIM>> geometry =
            new CartesianGridGeometry<NDIM>("geometry", app->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("hierarchy", geometry);
        Pointer<StandardTagAndInitialize<NDIM>> tagger =
            new StandardTagAndInitialize<NDIM>("tagger", ins, app->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer = new LoadBalancer<NDIM>("load_balancer");
        Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
            "gridding", app->getComponentDatabase("GriddingAlgorithm"), tagger, boxes, load_balancer);
        Pointer<CellVariable<NDIM, double>> ls = new CellVariable<NDIM, double>("ls");
        Pointer<CellVariable<NDIM, double>> H = new CellVariable<NDIM, double>("H");
        Pointer<CellVariable<NDIM, double>> T = new CellVariable<NDIM, double>("T");
        Pointer<CellVariable<NDIM, double>> lf = new CellVariable<NDIM, double>("lf");
        Pointer<CellVariable<NDIM, double>> extrap = new CellVariable<NDIM, double>("extrap");
        Pointer<CartGridFunction> initial =
            new muParserCartGridFunction("initial", app->getComponentDatabase("Initial"), geometry);
        for (auto var : { ls, H, T, lf, extrap })
        {
            adv->registerTransportedQuantity(var);
            adv->setAdvectionVelocity(var, ins->getAdvectionVelocityVariable());
            adv->setDiffusionCoefficient(var, 0.0);
            adv->setInitialConditions(var, initial);
        }
        Pointer<SideVariable<NDIM, double>> rho = new SideVariable<NDIM, double>("rho");
        ins->registerMassDensityVariable(rho);
        ins->registerMassDensityInitialConditions(initial);
        ins->initializePatchHierarchy(hierarchy, gridding);
        PIO::logOnlyNodeZero("output");
        plog << std::fixed << std::setprecision(12);
        const std::string contract = input->getString("contract");
        if (contract == "heaviside")
        {
            const int H_idx =
                VariableDatabase<NDIM>::getDatabase()->mapVariableAndContextToIndex(H, adv->getCurrentContext());
            check_heaviside(hierarchy, adv, ls, H_idx);
        }
        else if (contract == "capillary")
        {
            check_capillary(hierarchy, ins, adv, T, lf, extrap);
        }
        else
        {
            TBOX_ERROR("Unknown helper contract: " << contract << '\n');
        }
    }
    return 0;
}
