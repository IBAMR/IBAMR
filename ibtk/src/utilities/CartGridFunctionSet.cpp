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
CartGridFunctionSet::addFunction(SAMRAIPointer<CartGridFunction> fcn)
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
                                             SAMRAIPointer<VariableNd> var,
                                             SAMRAIPointer<PatchHierarchyNd> hierarchy,
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
    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
    const int cloned_data_idx = var_db->registerClonedPatchDataIndex(var, data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(cloned_data_idx);
    }
    SAMRAIPointer<HierarchyDataOpsRealNd<double> > hier_data_ops =
        HierarchyDataOpsManagerNd::getManager()->getOperationsDouble(var,
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
                                         SAMRAIPointer<VariableNd> var,
                                         SAMRAIPointer<PatchLevelNd> level,
                                         const double data_time,
                                         const bool initial_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level);
#endif
    SAMRAIPointer<CellVariableNd<double> > cc_var = var;
    SAMRAIPointer<EdgeVariableNd<double> > ec_var = var;
    SAMRAIPointer<FaceVariableNd<double> > fc_var = var;
    SAMRAIPointer<NodeVariableNd<double> > nc_var = var;
    SAMRAIPointer<SideVariableNd<double> > sc_var = var;
#if !defined(NDEBUG)
    TBOX_ASSERT(cc_var || ec_var || fc_var || nc_var || sc_var);
#endif
    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
    const int cloned_data_idx = var_db->registerClonedPatchDataIndex(var, data_idx);
    level->allocatePatchData(cloned_data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_fcns.empty());
#endif
    d_fcns[0]->setDataOnPatchLevel(data_idx, var, level, data_time, initial_time);
    for (unsigned int k = 1; k < d_fcns.size(); ++k)
    {
        d_fcns[k]->setDataOnPatchLevel(cloned_data_idx, var, level, data_time, initial_time);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = level->getPatch(p());
            if (cc_var)
            {
                SAMRAIPointer<CellDataNd<double> > data = patch->getPatchData(data_idx);
                SAMRAIPointer<CellDataNd<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchCellDataBasicOpsNd<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (ec_var)
            {
                SAMRAIPointer<EdgeDataNd<double> > data = patch->getPatchData(data_idx);
                SAMRAIPointer<EdgeDataNd<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchEdgeDataBasicOpsNd<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (fc_var)
            {
                SAMRAIPointer<FaceDataNd<double> > data = patch->getPatchData(data_idx);
                SAMRAIPointer<FaceDataNd<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchFaceDataBasicOpsNd<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (nc_var)
            {
                SAMRAIPointer<NodeDataNd<double> > data = patch->getPatchData(data_idx);
                SAMRAIPointer<NodeDataNd<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchNodeDataBasicOpsNd<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (sc_var)
            {
                SAMRAIPointer<SideDataNd<double> > data = patch->getPatchData(data_idx);
                SAMRAIPointer<SideDataNd<double> > cloned_data = patch->getPatchData(cloned_data_idx);
                PatchSideDataBasicOpsNd<double> patch_ops;
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
                                    SAMRAIPointer<VariableNd> var,
                                    SAMRAIPointer<PatchNd> patch,
                                    double data_time,
                                    bool initial_time,
                                    SAMRAIPointer<PatchLevelNd> patch_level)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(patch);
#endif
    SAMRAIPointer<CellVariableNd<double> > cc_var = var;
    SAMRAIPointer<EdgeVariableNd<double> > ec_var = var;
    SAMRAIPointer<FaceVariableNd<double> > fc_var = var;
    SAMRAIPointer<NodeVariableNd<double> > nc_var = var;
    SAMRAIPointer<SideVariableNd<double> > sc_var = var;
#if !defined(NDEBUG)
    TBOX_ASSERT(cc_var || ec_var || fc_var || nc_var || sc_var);
#endif
    SAMRAIPointer<PatchDataNd> data = patch->getPatchData(data_idx);
    SAMRAIPointer<PatchDataNd> cloned_data;
    if (cc_var)
    {
        SAMRAIPointer<CellDataNd<double> > p_data = data;
        cloned_data = new CellDataNd<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (ec_var)
    {
        SAMRAIPointer<EdgeDataNd<double> > p_data = data;
        cloned_data = new EdgeDataNd<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (fc_var)
    {
        SAMRAIPointer<FaceDataNd<double> > p_data = data;
        cloned_data = new FaceDataNd<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (nc_var)
    {
        SAMRAIPointer<NodeDataNd<double> > p_data = data;
        cloned_data = new NodeDataNd<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (sc_var)
    {
        SAMRAIPointer<SideDataNd<double> > p_data = data;
        cloned_data = new SideDataNd<double>(
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
            SAMRAIPointer<CellDataNd<double> > p_data = data;
            SAMRAIPointer<CellDataNd<double> > p_cloned_data = cloned_data;
            PatchCellDataBasicOpsNd<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (ec_var)
        {
            SAMRAIPointer<EdgeDataNd<double> > p_data = data;
            SAMRAIPointer<EdgeDataNd<double> > p_cloned_data = cloned_data;
            PatchEdgeDataBasicOpsNd<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (fc_var)
        {
            SAMRAIPointer<FaceDataNd<double> > p_data = data;
            SAMRAIPointer<FaceDataNd<double> > p_cloned_data = cloned_data;
            PatchFaceDataBasicOpsNd<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (nc_var)
        {
            SAMRAIPointer<NodeDataNd<double> > p_data = data;
            SAMRAIPointer<NodeDataNd<double> > p_cloned_data = cloned_data;
            PatchNodeDataBasicOpsNd<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (sc_var)
        {
            SAMRAIPointer<SideDataNd<double> > p_data = data;
            SAMRAIPointer<SideDataNd<double> > p_cloned_data = cloned_data;
            PatchSideDataBasicOpsNd<double> patch_ops;
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
