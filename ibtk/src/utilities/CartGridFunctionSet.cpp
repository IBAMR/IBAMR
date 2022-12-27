// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibtk/CartGridFunction.h"
#include "ibtk/CartGridFunctionSet.h"

#include "BasePatchLevel.h"
#include "CellData.h"
#include "CellVariable.h"
#include "EdgeData.h"
#include "EdgeVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "NodeData.h"
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchCellDataBasicOps.h"
#include "PatchData.h"
#include "PatchEdgeDataBasicOps.h"
#include "PatchFaceDataBasicOps.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchNodeDataBasicOps.h"
#include "PatchSideDataBasicOps.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

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
        if (fcn->isTimeDependent()) return true;
    }
    return false;
} // isTimeDependent

void
CartGridFunctionSet::setDataOnPatchHierarchy(const int data_idx,
                                             Pointer<Variable<NDIM> > var,
                                             Pointer<PatchHierarchy<NDIM> > hierarchy,
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
    Pointer<HierarchyDataOpsReal<NDIM, double> > hier_data_ops =
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
                                         Pointer<Variable<NDIM> > var,
                                         Pointer<PatchLevel<NDIM> > level,
                                         const double data_time,
                                         const bool initial_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level);
#endif
    Pointer<CellVariable<NDIM, double> > cc_var = var;
    Pointer<EdgeVariable<NDIM, double> > ec_var = var;
    Pointer<FaceVariable<NDIM, double> > fc_var = var;
    Pointer<NodeVariable<NDIM, double> > nc_var = var;
    Pointer<SideVariable<NDIM, double> > sc_var = var;
#if !defined(NDEBUG)
    TBOX_ASSERT(cc_var || ec_var || fc_var || nc_var || sc_var);
#endif
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
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            if (cc_var)
            {
                Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
                Pointer<CellData<NDIM, double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchCellDataBasicOps<NDIM, double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (ec_var)
            {
                Pointer<EdgeData<NDIM, double> > data = patch->getPatchData(data_idx);
                Pointer<EdgeData<NDIM, double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchEdgeDataBasicOps<NDIM, double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (fc_var)
            {
                Pointer<FaceData<NDIM, double> > data = patch->getPatchData(data_idx);
                Pointer<FaceData<NDIM, double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchFaceDataBasicOps<NDIM, double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (nc_var)
            {
                Pointer<NodeData<NDIM, double> > data = patch->getPatchData(data_idx);
                Pointer<NodeData<NDIM, double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchNodeDataBasicOps<NDIM, double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (sc_var)
            {
                Pointer<SideData<NDIM, double> > data = patch->getPatchData(data_idx);
                Pointer<SideData<NDIM, double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchSideDataBasicOps<NDIM, double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else
            {
                TBOX_ERROR(d_object_name << "::setDataOnPatchLevel():\n"
                                         << "  unsupported data centering.\n");
            }
        }
    }
    level->deallocatePatchData(cloned_data_idx);
    var_db->removePatchDataIndex(cloned_data_idx);
    return;
} // setDataOnPatchLevel

void
CartGridFunctionSet::setDataOnPatch(int data_idx,
                                    Pointer<Variable<NDIM> > var,
                                    Pointer<Patch<NDIM> > patch,
                                    double data_time,
                                    bool initial_time,
                                    Pointer<PatchLevel<NDIM> > patch_level)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(patch);
#endif
    Pointer<CellVariable<NDIM, double> > cc_var = var;
    Pointer<EdgeVariable<NDIM, double> > ec_var = var;
    Pointer<FaceVariable<NDIM, double> > fc_var = var;
    Pointer<NodeVariable<NDIM, double> > nc_var = var;
    Pointer<SideVariable<NDIM, double> > sc_var = var;
#if !defined(NDEBUG)
    TBOX_ASSERT(cc_var || ec_var || fc_var || nc_var || sc_var);
#endif
    Pointer<PatchData<NDIM> > data = patch->getPatchData(data_idx);
    Pointer<PatchData<NDIM> > cloned_data;
    if (cc_var)
    {
        Pointer<CellData<NDIM, double> > p_data = data;
        cloned_data = new CellData<NDIM, double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (ec_var)
    {
        Pointer<EdgeData<NDIM, double> > p_data = data;
        cloned_data = new EdgeData<NDIM, double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (fc_var)
    {
        Pointer<FaceData<NDIM, double> > p_data = data;
        cloned_data = new FaceData<NDIM, double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (nc_var)
    {
        Pointer<NodeData<NDIM, double> > p_data = data;
        cloned_data = new NodeData<NDIM, double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (sc_var)
    {
        Pointer<SideData<NDIM, double> > p_data = data;
        cloned_data = new SideData<NDIM, double>(
            p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth(), p_data->getDirectionVector());
    }
    else
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                                 << "  unsupported data centering.\n");
    }
    cloned_data->setTime(data->getTime());
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_fcns.empty());
#endif
    d_fcns[0]->setDataOnPatch(data_idx, var, patch, data_time, initial_time, patch_level);
    cloned_data->copy(*data);
    // NOTE: We operate on data_idx instead of cloned_data_idx here because it
    // is not straightforward to add a cloned data index to a single patch.
    for (unsigned int k = 1; k < d_fcns.size(); ++k)
    {
        d_fcns[k]->setDataOnPatch(data_idx, var, patch, data_time, initial_time, patch_level);
        if (cc_var)
        {
            Pointer<CellData<NDIM, double> > p_data = data;
            Pointer<CellData<NDIM, double> > p_cloned_data = cloned_data;
            PatchCellDataBasicOps<NDIM, double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (ec_var)
        {
            Pointer<EdgeData<NDIM, double> > p_data = data;
            Pointer<EdgeData<NDIM, double> > p_cloned_data = cloned_data;
            PatchEdgeDataBasicOps<NDIM, double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (fc_var)
        {
            Pointer<FaceData<NDIM, double> > p_data = data;
            Pointer<FaceData<NDIM, double> > p_cloned_data = cloned_data;
            PatchFaceDataBasicOps<NDIM, double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (nc_var)
        {
            Pointer<NodeData<NDIM, double> > p_data = data;
            Pointer<NodeData<NDIM, double> > p_cloned_data = cloned_data;
            PatchNodeDataBasicOps<NDIM, double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (sc_var)
        {
            Pointer<SideData<NDIM, double> > p_data = data;
            Pointer<SideData<NDIM, double> > p_cloned_data = cloned_data;
            PatchSideDataBasicOps<NDIM, double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else
        {
            TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                                     << "  unsupported data centering.\n");
        }
    }
    data->copy(*cloned_data);
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
