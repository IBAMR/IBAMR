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

#include <ibtk/AppInitializer.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Logger.h>

#include <CartesianGridGeometry.h>
#include <PatchHierarchy.h>
#include <SideGeometry.h>
#include <SideVariable.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <utility>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

class ObservedGeometry : public CartesianGridGeometry<NDIM>
{
public:
    using CartesianGridGeometry<NDIM>::CartesianGridGeometry;

    void setGeometryDataOnPatch(Patch<NDIM>& patch,
                                const IntVector<NDIM>& ratio,
                                const tbox::Array<tbox::Array<bool>>& regular,
                                const tbox::Array<tbox::Array<bool>>& periodic) const override
    {
        ++d_created;
        CartesianGridGeometry<NDIM>::setGeometryDataOnPatch(patch, ratio, regular, periodic);
    }

    int getCreatedCount() const
    {
        return d_created;
    }

private:
    mutable int d_created = 0;
};

class PolicyTestAppender : public TestAppender
{
public:
    explicit PolicyTestAppender(std::function<bool()> unchanged) : d_unchanged(std::move(unchanged))
    {
    }
    void logMessage(const std::string& message, const std::string& file, const int line) override
    {
        plog << "Patch data unchanged before rejection: " << (d_unchanged() ? "true" : "false") << std::endl;
        TestAppender::logMessage(message, file, line);
    }

private:
    std::function<bool()> d_unchanged;
};

std::string
capture_boundary_boxes(Pointer<PatchHierarchy<NDIM>> hierarchy, const bool require_wall_boxes)
{
    std::ostringstream boxes;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<PatchGeometry<NDIM>> geometry = level->getPatch(p())->getPatchGeometry();
            // Every patch in this nonperiodic fixture meets a physical wall.
            if (require_wall_boxes && geometry->getCodimensionBoundaries(1).size() == 0)
            {
                TBOX_ERROR("Deferred boundary boxes were not constructed on every required level.\n");
            }
            for (int codim = 1; codim <= NDIM; ++codim)
            {
                const auto& boundaries = geometry->getCodimensionBoundaries(codim);
                for (int k = 0; k < boundaries.size(); ++k)
                {
                    boxes << ln << ' ' << p() << ' ' << codim << ' ' << boundaries[k].getLocationIndex() << ' '
                          << boundaries[k].getBox() << '\n';
                }
            }
        }
    }
    return boxes.str();
}

double
exact_value(Pointer<Patch<NDIM>> patch,
            const SideIndex<NDIM>& index,
            const int axis,
            const bool constant,
            const bool quadratic)
{
    if (constant)
    {
        return 3.0 + axis;
    }
    Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    double value = 3.0 + axis;
    for (int d = 0; d < NDIM; ++d)
    {
        const double x = geometry->getXLower()[d] +
                         geometry->getDx()[d] * (index(d) - patch->getBox().lower()(d) + (d == axis ? 0.0 : 0.5));
        value += (d + 1) * (quadratic ? x * x : x);
    }
    return value;
}

void
seed(Pointer<PatchHierarchy<NDIM>> hierarchy, const int data_idx, const bool constant, const bool quadratic)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(data_idx))
        {
            level->allocatePatchData(data_idx);
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> data = patch->getPatchData(data_idx);
            data->fillAll(std::numeric_limits<double>::quiet_NaN());
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> i(patch->getBox(), axis); i; i++)
                {
                    (*data)(i()) = exact_value(patch, i(), axis, constant, quadratic);
                }
            }
        }
    }
}

bool
unchanged(Pointer<PatchHierarchy<NDIM>> hierarchy, const int data_idx, const bool constant, const bool quadratic)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        // Only the transaction data were allocated before initialization.
        // In particular, a rejected internal registration must not allocate
        // the coarse-fine helper's scratch data before the width check.
        const int n_components = hierarchy->getPatchDescriptor()->getMaxNumberRegisteredComponents();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            for (int idx = 0; idx < n_components; ++idx)
            {
                if (idx != data_idx && patch->checkAllocated(idx))
                {
                    return false;
                }
            }
            Pointer<SideData<NDIM, double>> data = patch->getPatchData(data_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> interior = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
                for (SideIterator<NDIM> i(data->getGhostBox(), axis); i; i++)
                {
                    const double value = (*data)(i());
                    if (interior.contains(i()))
                    {
                        if (value != exact_value(patch, i(), axis, constant, quadratic))
                        {
                            return false;
                        }
                    }
                    else if (!std::isnan(value))
                    {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

void
check(Pointer<PatchHierarchy<NDIM>> hierarchy,
      const int data_idx,
      const bool constant,
      const bool quadratic,
      const int coarsest)
{
    int bad = 0;
    double error = 0.0;
    for (int ln = coarsest; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> data = patch->getPatchData(data_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> i(data->getGhostBox(), axis); i; i++)
                {
                    const double value = (*data)(i());
                    if (!std::isfinite(value))
                    {
                        ++bad;
                    }
                    else
                    {
                        // Check quadratic data only on faces with domain cells on
                        // both sides. Linear physical boundary extrapolation need
                        // not reproduce quadratic data on boundary faces.
                        hier::Index<NDIM> lower_cell = i(), upper_cell = i();
                        --lower_cell(axis);
                        bool lower_inside = false, upper_inside = false;
                        const BoxArray<NDIM>& domain = level->getPhysicalDomain();
                        for (int k = 0; k < domain.size(); ++k)
                        {
                            lower_inside = lower_inside || domain[k].contains(lower_cell);
                            upper_inside = upper_inside || domain[k].contains(upper_cell);
                        }
                        if (quadratic && !(lower_inside && upper_inside))
                        {
                            continue;
                        }
                        const double expected = exact_value(patch, i(), axis, constant, quadratic);
                        const double local_error = std::abs(value - expected);
                        error = std::max(error, local_error);
                    }
                }
            }
        }
    }
    bad = IBTK_MPI::sumReduction(bad);
    error = IBTK_MPI::maxReduction(error);
    plog << "Nonfinite entries: " << bad << "; maximum error: " << error << std::endl;
    if (bad || error > 1.0e-10)
    {
        TBOX_ERROR("Ghost-width policy test failed field verification.\n");
    }
}

int
register_side_variable(const std::string& name, const IntVector<NDIM>& width)
{
    auto* db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double>> variable = new SideVariable<NDIM, double>(name);
    return db->registerVariableAndContext(variable, db->getContext("ghost_width_policy"), width);
}

int
report_missing_rejection()
{
    plog << "Expected width-policy rejection was not raised." << std::endl;
    Logger::getInstance()->setAbortAppender(new TestAppender());
    return 0;
}

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    {
        TimerManager::createManager(nullptr);
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
        Pointer<Database> input = app_initializer->getInputDatabase();
        const std::string mode = input->getString("test");
        plog << "Case: " << mode << std::endl;

        const bool deferred = mode == "deferred-before" || mode == "deferred-after" || mode == "deferred-range";
        const bool partial_range = mode == "range" || mode == "deferred-range";
        const bool periodic = mode == "periodic" || mode == "periodic-late";
        const bool use_cf = mode == "cf-constant" || mode == "cf-late-scratch";
        const bool constant = periodic || use_cf;
        const bool quadratic = mode == "l-shape";
        const bool multilevel = mode != "single" && mode != "late-init" && mode != "deferred-after" &&
                                mode != "cf-late-scratch" && !quadratic;
        const bool anisotropic = mode == "anisotropic-valid" || mode == "anisotropic-reject";
        const bool reserved = mode == "shrink-restore" || mode == "shrink-fill-restore" || mode == "temporary-shrink" ||
                              mode == "restore-unobserved" || mode == "retained-reservation" || anisotropic;
        const int initial_width = mode == "cf-late-scratch" ? 0 :
                                  mode == "early" || mode == "single" || partial_range || mode == "coarsen" ||
                                          mode == "fresh-geometry" || periodic || use_cf || quadratic ?
                                                              2 :
                                                              1;
        auto* db = VariableDatabase<NDIM>::getDatabase();
        int data_idx = register_side_variable("field", IntVector<NDIM>(initial_width));
        int reserve_idx = -1;
        if (reserved)
        {
            IntVector<NDIM> width(anisotropic ? 1 : 2);
            width(0) = 2;
            reserve_idx = register_side_variable("unused_reservation", width);
        }

        BoxArray<NDIM> domain(quadratic ? 2 : 1), coarse(2), fine(2);
        domain[0] = Box<NDIM>(hier::Index<NDIM>(0), hier::Index<NDIM>(15));
        coarse[0] = coarse[1] = domain[0];
        coarse[0].upper()(0) = 7;
        coarse[1].lower()(0) = 8;
        if (quadratic)
        {
            coarse[1].upper()(1) = 7;
            domain = coarse;
        }
        fine[0] = fine[1] = Box<NDIM>(hier::Index<NDIM>(0), hier::Index<NDIM>(31));
        fine[0].lower()(0) = 8;
        fine[0].upper()(0) = 15;
        fine[1].lower()(0) = 16;
        fine[1].upper()(0) = 23;
        fine[0].upper()(1) = fine[1].upper()(1) = 15;
        double lo[NDIM], hi[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            lo[d] = 0.0;
            hi[d] = 16.0;
        }
        Pointer<ObservedGeometry> geometry = new ObservedGeometry("geometry", lo, hi, domain, false);
        if (periodic)
        {
            geometry->initializePeriodicShift(IntVector<NDIM>(1));
        }
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("hierarchy", geometry, false);
        using Transaction = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        if (mode == "empty")
        {
            HierarchyGhostCellInterpolation empty_op;
            Transaction empty_transaction(data_idx, "CONSERVATIVE_LINEAR_REFINE", false, "NONE", "LINEAR");
            empty_op.initializeOperatorState(empty_transaction, hierarchy);
            empty_op.resetTransactionComponent(empty_transaction);
            empty_op.fillData(0.0);
            register_side_variable("before_first_level", IntVector<NDIM>(2));
        }
        ProcessorMapping mapping(2);
        mapping.setProcessorAssignment(0, 0);
        mapping.setProcessorAssignment(1, IBTK_MPI::getNodes() > 1 ? 1 : 0);
        hierarchy->makeNewPatchLevel(0, IntVector<NDIM>(1), coarse, mapping, deferred);
        if (multilevel)
        {
            hierarchy->makeNewPatchLevel(1, IntVector<NDIM>(2), fine, mapping, deferred);
        }
        if (mode == "deferred-before" || mode == "late-init")
        {
            data_idx = register_side_variable("before_init", IntVector<NDIM>(2));
        }
        const std::string original_boundary_boxes = capture_boundary_boxes(hierarchy, !deferred && !periodic);
        if (mode == "shrink-restore" || mode == "temporary-shrink")
        {
            db->removePatchDataIndex(reserve_idx);
        }
        if (mode == "temporary-shrink")
        {
            // Native construction can lower the shared limit even before HGC
            // observes the narrower descriptor. Destroying that level cannot
            // restore the limit used by the older, surviving boundary boxes.
            {
                Pointer<PatchLevel<NDIM>> temporary = new PatchLevel<NDIM>(
                    coarse, mapping, IntVector<NDIM>(1), geometry, hierarchy->getPatchDescriptor());
            }
            register_side_variable("replacement", IntVector<NDIM>(2));
        }
        seed(hierarchy, data_idx, constant, quadratic);
        Pointer<Logger::Appender> appender =
            new PolicyTestAppender([&]() { return unchanged(hierarchy, data_idx, constant, quadratic); });
        Logger::getInstance()->setAbortAppender(appender);

        Transaction transaction(data_idx,
                                mode == "cf-late-scratch" ? "NONE" : "CONSERVATIVE_LINEAR_REFINE",
                                use_cf,
                                mode == "coarsen" ? "CONSERVATIVE_COARSEN" : "NONE",
                                "LINEAR");
        const int coarsest = partial_range ? 1 : 0;
        const int finest = hierarchy->getFinestLevelNumber();
        const int before = geometry->getCreatedCount();
        HierarchyGhostCellInterpolation op;
        plog << "Checking initialization." << std::endl;
        op.initializeOperatorState(transaction, hierarchy, coarsest, finest);
        if (mode == "late-init" || mode == "temporary-shrink" || mode == "cf-late-scratch")
        {
            return report_missing_rejection();
        }
        const std::string initialized_boundary_boxes = capture_boundary_boxes(hierarchy, !periodic);
        if (!deferred && original_boundary_boxes != initialized_boundary_boxes)
        {
            TBOX_ERROR("Initialization changed existing boundary boxes.\n");
        }
        const int created = IBTK_MPI::sumReduction(geometry->getCreatedCount() - before);
        if (multilevel && created == 0)
        {
            TBOX_ERROR("Expected temporary coarse-level geometry.\n");
        }

        if (mode == "shrink-restore")
        {
            op.deallocateOperatorState();
            register_side_variable("replacement", IntVector<NDIM>(2));
            plog << "Checking reinitialization after temporary levels were destroyed." << std::endl;
            op.initializeOperatorState(transaction, hierarchy, coarsest, finest);
            return report_missing_rejection();
        }
        if (mode == "shrink-fill-restore" || mode == "restore-unobserved")
        {
            db->removePatchDataIndex(reserve_idx);
        }
        if (mode == "late-fill" || mode == "late-reset" || mode == "deferred-after" || mode == "periodic-late" ||
            mode == "restore-unobserved" || mode == "retained-reservation")
        {
            register_side_variable("later", IntVector<NDIM>(periodic ? 3 : 2));
        }
        if (mode == "retained-reservation")
        {
            register_side_variable("later_smaller", IntVector<NDIM>(1));
        }
        if (anisotropic)
        {
            IntVector<NDIM> width(1);
            width(mode == "anisotropic-valid" ? 0 : NDIM - 1) = 2;
            register_side_variable("later_anisotropic", width);
        }
        if (mode == "late-reset")
        {
            plog << "Checking reset." << std::endl;
            op.resetTransactionComponent(transaction);
            return report_missing_rejection();
        }
        plog << "Checking fill." << std::endl;
        op.fillData(0.0);
        if (mode == "late-fill" || mode == "deferred-after" || mode == "periodic-late" || mode == "anisotropic-reject")
        {
            return report_missing_rejection();
        }
        check(hierarchy, data_idx, constant, quadratic, coarsest);
        if (quadratic)
        {
            // A subsequent exchange can mask overwrites inside the domain.
            // Also check the physical boundary operator on already filled data.
            CartExtrapPhysBdryOp extrapolation(data_idx, "LINEAR");
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                extrapolation.setPhysicalBoundaryConditions(*level->getPatch(p()), 0.0, IntVector<NDIM>(2));
            }
            check(hierarchy, data_idx, constant, quadratic, coarsest);
        }
        if (mode == "shrink-fill-restore")
        {
            seed(hierarchy, data_idx, constant, quadratic);
            register_side_variable("replacement", IntVector<NDIM>(2));
            plog << "Checking fill after the limit decreased." << std::endl;
            op.fillData(0.0);
            return report_missing_rejection();
        }

        op.resetTransactionComponent(transaction);
        seed(hierarchy, data_idx, constant, quadratic);
        op.fillData(0.0);
        check(hierarchy, data_idx, constant, quadratic, coarsest);
        op.reinitializeOperatorState(hierarchy);
        seed(hierarchy, data_idx, constant, quadratic);
        op.fillData(0.0);
        check(hierarchy, data_idx, constant, quadratic, 0);
        HierarchyGhostCellInterpolation other;
        other.initializeOperatorState(transaction, hierarchy);
        op.deallocateOperatorState();
        seed(hierarchy, data_idx, constant, quadratic);
        other.fillData(0.0);
        check(hierarchy, data_idx, constant, quadratic, 0);

        if (mode == "fresh-geometry")
        {
            register_side_variable("wider_on_new_geometry", IntVector<NDIM>(3));
            geometry = new ObservedGeometry("new_geometry", lo, hi, domain, false);
            hierarchy = new PatchHierarchy<NDIM>("new_hierarchy", geometry, false);
            hierarchy->makeNewPatchLevel(0, IntVector<NDIM>(1), coarse, mapping);
            hierarchy->makeNewPatchLevel(1, IntVector<NDIM>(2), fine, mapping);
            seed(hierarchy, data_idx, constant, quadratic);
            op.initializeOperatorState(transaction, hierarchy);
            op.fillData(0.0);
            check(hierarchy, data_idx, constant, quadratic, 0);
        }
        plog << "Completed." << std::endl;
        Logger::getInstance()->setAbortAppender(new TestAppender());
    }
}
