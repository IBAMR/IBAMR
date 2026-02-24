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
#include "ibtk/samrai_compatibility_names.h"

#include "MultiblockDataTranslator.h"
#include "SAMRAIBasePatchLevel.h"
#include "SAMRAICellData.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIEdgeData.h"
#include "SAMRAIEdgeVariable.h"
#include "SAMRAIFaceData.h"
#include "SAMRAIFaceVariable.h"
#include "SAMRAIHierarchyDataOpsManager.h"
#include "SAMRAIHierarchyDataOpsReal.h"
#include "SAMRAIIntVector.h"
#include "SAMRAINodeData.h"
#include "SAMRAINodeVariable.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchCellDataBasicOps.h"
#include "SAMRAIPatchData.h"
#include "SAMRAIPatchEdgeDataBasicOps.h"
#include "SAMRAIPatchFaceDataBasicOps.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPatchNodeDataBasicOps.h"
#include "SAMRAIPatchSideDataBasicOps.h"
#include "SAMRAIPointer.h"
#include "SAMRAISideData.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableDatabase.h"

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
                                             SAMRAIPointer<SAMRAIVariable> var,
                                             SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
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
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    const int cloned_data_idx = var_db->registerClonedPatchDataIndex(var, data_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(cloned_data_idx);
    }
    SAMRAIPointer<SAMRAIHierarchyDataOpsReal<double>> hier_data_ops =
        SAMRAIHierarchyDataOpsManager::getManager()->getOperationsDouble(var,
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
                                         SAMRAIPointer<SAMRAIVariable> var,
                                         SAMRAIPointer<SAMRAIPatchLevel> level,
                                         const double data_time,
                                         const bool initial_time)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level);
#endif
    SAMRAIPointer<SAMRAICellVariable<double>> cc_var = var;
    SAMRAIPointer<SAMRAIEdgeVariable<double>> ec_var = var;
    SAMRAIPointer<SAMRAIFaceVariable<double>> fc_var = var;
    SAMRAIPointer<SAMRAINodeVariable<double>> nc_var = var;
    SAMRAIPointer<SAMRAISideVariable<double>> sc_var = var;
#if !defined(NDEBUG)
    TBOX_ASSERT(cc_var || ec_var || fc_var || nc_var || sc_var);
#endif
    SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
    const int cloned_data_idx = var_db->registerClonedPatchDataIndex(var, data_idx);
    level->allocatePatchData(cloned_data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_fcns.empty());
#endif
    d_fcns[0]->setDataOnPatchLevel(data_idx, var, level, data_time, initial_time);
    for (unsigned int k = 1; k < d_fcns.size(); ++k)
    {
        d_fcns[k]->setDataOnPatchLevel(cloned_data_idx, var, level, data_time, initial_time);
        for (SAMRAIPatchLevel::Iterator p(level); p; p++)
        {
            SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
            if (cc_var)
            {
                SAMRAIPointer<SAMRAICellData<double>> data = patch->getPatchData(data_idx);
                SAMRAIPointer<SAMRAICellData<double>> cloned_data = patch->getPatchData(cloned_data_idx);
                SAMRAIPatchCellDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (ec_var)
            {
                SAMRAIPointer<SAMRAIEdgeData<double>> data = patch->getPatchData(data_idx);
                SAMRAIPointer<SAMRAIEdgeData<double>> cloned_data = patch->getPatchData(cloned_data_idx);
                SAMRAIPatchEdgeDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (fc_var)
            {
                SAMRAIPointer<SAMRAIFaceData<double>> data = patch->getPatchData(data_idx);
                SAMRAIPointer<SAMRAIFaceData<double>> cloned_data = patch->getPatchData(cloned_data_idx);
                SAMRAIPatchFaceDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (nc_var)
            {
                SAMRAIPointer<SAMRAINodeData<double>> data = patch->getPatchData(data_idx);
                SAMRAIPointer<SAMRAINodeData<double>> cloned_data = patch->getPatchData(cloned_data_idx);
                SAMRAIPatchNodeDataBasicOps<double> patch_ops;
                patch_ops.add(data, data, cloned_data, patch->getBox());
            }
            else if (sc_var)
            {
                SAMRAIPointer<SAMRAISideData<double>> data = patch->getPatchData(data_idx);
                SAMRAIPointer<SAMRAISideData<double>> cloned_data = patch->getPatchData(cloned_data_idx);
                SAMRAIPatchSideDataBasicOps<double> patch_ops;
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
                                    SAMRAIPointer<SAMRAIVariable> var,
                                    SAMRAIPointer<SAMRAIPatch> patch,
                                    double data_time,
                                    bool initial_time,
                                    SAMRAIPointer<SAMRAIPatchLevel> patch_level)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(patch);
#endif
    SAMRAIPointer<SAMRAICellVariable<double>> cc_var = var;
    SAMRAIPointer<SAMRAIEdgeVariable<double>> ec_var = var;
    SAMRAIPointer<SAMRAIFaceVariable<double>> fc_var = var;
    SAMRAIPointer<SAMRAINodeVariable<double>> nc_var = var;
    SAMRAIPointer<SAMRAISideVariable<double>> sc_var = var;
#if !defined(NDEBUG)
    TBOX_ASSERT(cc_var || ec_var || fc_var || nc_var || sc_var);
#endif
    SAMRAIPointer<SAMRAIPatchData> data = patch->getPatchData(data_idx);
    SAMRAIPointer<SAMRAIPatchData> cloned_data;
    if (cc_var)
    {
        SAMRAIPointer<SAMRAICellData<double>> p_data = data;
        cloned_data = new SAMRAICellData<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (ec_var)
    {
        SAMRAIPointer<SAMRAIEdgeData<double>> p_data = data;
        cloned_data = new SAMRAIEdgeData<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (fc_var)
    {
        SAMRAIPointer<SAMRAIFaceData<double>> p_data = data;
        cloned_data = new SAMRAIFaceData<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (nc_var)
    {
        SAMRAIPointer<SAMRAINodeData<double>> p_data = data;
        cloned_data = new SAMRAINodeData<double>(p_data->getBox(), p_data->getDepth(), p_data->getGhostCellWidth());
    }
    else if (sc_var)
    {
        SAMRAIPointer<SAMRAISideData<double>> p_data = data;
        cloned_data = new SAMRAISideData<double>(
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
            SAMRAIPointer<SAMRAICellData<double>> p_data = data;
            SAMRAIPointer<SAMRAICellData<double>> p_cloned_data = cloned_data;
            SAMRAIPatchCellDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (ec_var)
        {
            SAMRAIPointer<SAMRAIEdgeData<double>> p_data = data;
            SAMRAIPointer<SAMRAIEdgeData<double>> p_cloned_data = cloned_data;
            SAMRAIPatchEdgeDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (fc_var)
        {
            SAMRAIPointer<SAMRAIFaceData<double>> p_data = data;
            SAMRAIPointer<SAMRAIFaceData<double>> p_cloned_data = cloned_data;
            SAMRAIPatchFaceDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (nc_var)
        {
            SAMRAIPointer<SAMRAINodeData<double>> p_data = data;
            SAMRAIPointer<SAMRAINodeData<double>> p_cloned_data = cloned_data;
            SAMRAIPatchNodeDataBasicOps<double> patch_ops;
            patch_ops.add(p_cloned_data, p_cloned_data, p_data, patch->getBox());
        }
        else if (sc_var)
        {
            SAMRAIPointer<SAMRAISideData<double>> p_data = data;
            SAMRAIPointer<SAMRAISideData<double>> p_cloned_data = cloned_data;
            SAMRAIPatchSideDataBasicOps<double> patch_ops;
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
