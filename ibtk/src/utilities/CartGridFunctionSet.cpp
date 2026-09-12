// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/CartGridFunction.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/CartesianCentering.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <BasePatchLevel.h>
#include <HierarchyDataOpsManager.h>
#include <HierarchyDataOpsReal.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <Variable.h>
#include <VariableDatabase.h>

#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartGridFunctionSet::CartGridFunctionSet(std::string object_name) : CartGridFunction(std::move(object_name))
{
    // intentionally blank
    return;
} // CartGridFunctionSet

void
CartGridFunctionSet::addFunction(Pointer<CartGridFunction> fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(fcn);
#endif
    d_fcns.push_back(fcn);
    return;
} // addFunction

bool
CartGridFunctionSet::isTimeDependent() const
{
    for (const auto& fcn : d_fcns)
    {
        if (fcn->isTimeDependent())
        {
            return true;
        }
    }
    return false;
} // isTimeDependent

void
CartGridFunctionSet::setDataOnPatchHierarchy(const int data_idx,
                                             Pointer<Variable<NDIM>> var,
                                             Pointer<PatchHierarchy<NDIM>> hierarchy,
                                             const double data_time,
                                             const bool initial_time,
                                             const int coarsest_ln_in,
                                             const int finest_ln_in)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
#endif
    const int coarsest_ln = (coarsest_ln_in == invalid_level_number ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln_in);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int cloned_data_idx = var_db->registerClonedPatchDataIndex(var, data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(cloned_data_idx);
    }
    Pointer<HierarchyDataOpsReal<NDIM, double>> hier_data_ops =
        HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(var,
                                                                         hierarchy,
                                                                         /* get_unique */ true);
    if (!hier_data_ops)
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatchHierarchy():\n"
                                 << "  unsupported data centering.\n");
    }
    hier_data_ops->resetLevels(coarsest_ln, finest_ln);
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_fcns.empty());
#endif
    d_fcns[0]->setDataOnPatchHierarchy(data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);
    for (unsigned int k = 1; k < d_fcns.size(); ++k)
    {
        d_fcns[k]->setDataOnPatchHierarchy(
            cloned_data_idx, var, hierarchy, data_time, initial_time, coarsest_ln_in, finest_ln_in);
        hier_data_ops->add(data_idx, data_idx, cloned_data_idx);
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(cloned_data_idx);
    }
    var_db->removePatchDataIndex(cloned_data_idx);
    return;
} // setDataOnPatchHierarchy

void
CartGridFunctionSet::setDataOnPatchLevel(const int data_idx,
                                         Pointer<Variable<NDIM>> var,
                                         Pointer<PatchLevel<NDIM>> level,
                                         const double data_time,
                                         const bool initial_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level);
#endif
    const auto add_functions = [&]<Centering Layout>()
    {
        using Data = typename Layout::template Data<double>;
        typename Layout::template PatchOps<double> patch_ops;
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int cloned_data_idx = var_db->registerClonedPatchDataIndex(var, data_idx);
        level->allocatePatchData(cloned_data_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(!d_fcns.empty());
#endif
        d_fcns[0]->setDataOnPatchLevel(data_idx, var, level, data_time, initial_time);
        for (unsigned int k = 1; k < d_fcns.size(); ++k)
        {
            d_fcns[k]->setDataOnPatchLevel(cloned_data_idx, var, level, data_time, initial_time);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<Data> data = patch->getPatchData(data_idx);
                Pointer<Data> cloned_data = patch->getPatchData(cloned_data_idx);
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
        }
        level->deallocatePatchData(cloned_data_idx);
        var_db->removePatchDataIndex(cloned_data_idx);
    };
    dispatch_data_centering(get_data_centering<double>(*var->getPatchDataFactory()), add_functions);
} // setDataOnPatchLevel

void
CartGridFunctionSet::setDataOnPatch(const int data_idx,
                                    Pointer<Variable<NDIM>> var,
                                    Pointer<Patch<NDIM>> patch,
                                    const double data_time,
                                    const bool initial_time,
                                    Pointer<PatchLevel<NDIM>> patch_level)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(patch);
#endif
    dispatch_data_centering(
        get_data_centering<double>(*var->getPatchDataFactory()),
        [&]<Centering Layout>()
        {
            using Data = typename Layout::template Data<double>;
            typename Layout::template PatchOps<double> patch_ops;
            Pointer<Data> data = patch->getPatchData(data_idx);
            Pointer<Data> cloned_data;
            if constexpr (requires { data->getDirectionVector(); })
            {
                cloned_data =
                    new Data(data->getBox(), data->getDepth(), data->getGhostCellWidth(), data->getDirectionVector());
            }
            else
            {
                cloned_data = new Data(data->getBox(), data->getDepth(), data->getGhostCellWidth());
            }
            cloned_data->setTime(data->getTime());
#if !defined(NDEBUG)
            TBOX_ASSERT(!d_fcns.empty());
#endif
            d_fcns[0]->setDataOnPatch(data_idx, var, patch, data_time, initial_time, patch_level);
            cloned_data->copy(*data);
            // A cloned index cannot be registered for a single patch, so accumulate
            // separately while each function writes to the original index.
            for (unsigned int k = 1; k < d_fcns.size(); ++k)
            {
                d_fcns[k]->setDataOnPatch(data_idx, var, patch, data_time, initial_time, patch_level);
                patch_ops.add(cloned_data, cloned_data, data, patch->getBox());
            }
            data->copy(*cloned_data);
        });
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
