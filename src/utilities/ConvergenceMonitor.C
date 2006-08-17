// Filename: ConvergenceMonitor.C
// Last modified: <17.Aug.2006 16:03:00 boyce@bigboy.nyconnect.com>
// Created on 19 Jun 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ConvergenceMonitor.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

#include <Box.h>
#include <BoxArray.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <HierarchyDataOpsManager.h>
#include <HierarchyDataOpsReal.h>
#include <IntVector.h>
#include <Patch.h>
#include <VariableDatabase.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

ConvergenceMonitor::ConvergenceMonitor(
    const std::string& object_name)
{
    d_object_name = object_name;

    // Initialize Variables and contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    
    d_context = var_db->getContext(d_object_name+"::CONTEXT");
    d_wgt_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::wgt",1);
    d_wgt_idx = var_db->
        registerVariableAndContext(d_wgt_var, d_context, SAMRAI::hier::IntVector<NDIM>(0));
    
    return;
}// ConvergenceMonitor

ConvergenceMonitor::~ConvergenceMonitor()
{
    // intentionally blank
    return;
}// ~ConvergenceMonitor

void
ConvergenceMonitor::registerMonitoredVariableAndContext(
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> var_ctx,
    SAMRAI::tbox::Pointer<SetDataStrategy> exact_soln_setter)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!var.isNull());
    assert(!var_ctx.isNull());
    assert(!exact_soln_setter.isNull());
#endif

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    const int var_idx = var_db->mapVariableAndContextToIndex(var, var_ctx);

    d_monitored_indices.insert(var_idx);
    d_monitored_vars[var_idx] = var;
    d_monitored_var_ctxs[var_idx] = var_ctx;
    d_exact_soln_setters[var_idx] = exact_soln_setter;
    return;
}// registerMonitoredVariableAndContext

void
ConvergenceMonitor::monitorConvergence(
    const double data_time)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_data_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();
    
    for (set<int>::const_iterator it = d_monitored_indices.begin();
         it != d_monitored_indices.end(); ++it)
    {
        const int var_idx = *it;
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var = d_monitored_vars[var_idx];
        SAMRAI::tbox::Pointer<SetDataStrategy>& exact_soln_setter =
            d_exact_soln_setters[var_idx];
        
        const int cloned_idx = var_db->
            registerClonedPatchDataIndex(var,var_idx);
        
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->
                allocatePatchData(cloned_idx,data_time);
        }

        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM,double> > hier_data_ops =
            hier_data_ops_manager->getOperationsDouble(var, d_hierarchy);
        hier_data_ops->resetLevels(d_coarsest_ln, d_finest_ln);

        // Compute the exact solution.
        exact_soln_setter->setDataOnPatchHierarchy(
            cloned_idx, var, d_hierarchy, data_time);
        
        hier_data_ops->subtract(cloned_idx, var_idx, cloned_idx);
        
        // Compute error norms.
        //
        // For cell centered variables, we weight each cell by its
        // length/area/volume when computing error norms.
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = var;
        
        SAMRAI::tbox::pout << "Error in " << var->getName()
                           << " at time " << data_time << ":\n";
        
        if (cc_var.isNull())
        {
            SAMRAI::tbox::pout << "  L1-norm:  " << hier_data_ops->L1Norm(cloned_idx)  << "\n";
            SAMRAI::tbox::pout << "  L2-norm:  " << hier_data_ops->L2Norm(cloned_idx)  << "\n";
            SAMRAI::tbox::pout << "  max-norm: " << hier_data_ops->maxNorm(cloned_idx) << "\n";
        }
        else
        {
            SAMRAI::tbox::pout << "  L1-norm:  " << hier_data_ops->L1Norm(cloned_idx,d_wgt_idx)  << "\n";
            SAMRAI::tbox::pout << "  L2-norm:  " << hier_data_ops->L2Norm(cloned_idx,d_wgt_idx)  << "\n";
            SAMRAI::tbox::pout << "  max-norm: " << hier_data_ops->maxNorm(cloned_idx,d_wgt_idx) << "\n";
        }
        
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(cloned_idx);
        }
        
        var_db->removePatchDataIndex(cloned_idx);
    }
    return;
}// monitorConvergence

void
ConvergenceMonitor::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0) 
           && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!(old_level.isNull())) {
        assert(level_number == old_level->getLevelNumber());
    }
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    // Allocate all patch data managed by the monitor.
    if (allocate_data)
    {
        hierarchy->getPatchLevel(level_number)->allocatePatchData(d_wgt_idx);
    }
    return;
}// initializeLevelData

void
ConvergenceMonitor::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
           && (coarsest_level <= finest_level) 
           && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln0 = 0; ln0 <= finest_level; ++ln0)
    {
        assert(!(hierarchy->getPatchLevel(ln0)).isNull());
    }
#endif
    // Reset the Hierarchy.
    d_hierarchy = hierarchy;
    d_coarsest_ln = 0;
    d_finest_ln = d_hierarchy->getFinestLevelNumber();
    
    // Each cell's weight is set to its cell volume, unless the cell
    // is refined on a finer level, in which case the weight is set to
    // zero.  This insures that no part of the physical domain is
    // counted twice when discrete norms and integrals are calculated
    // on the entire hierarchy.
    const int finest_hier_level = d_finest_ln;
    
    for (int ln = SAMRAI::tbox::Utilities::imax(coarsest_level-1,0);
         ln <= finest_hier_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        SAMRAI::hier::BoxArray<NDIM> refined_region_boxes;
        
        if (ln < finest_hier_level)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > next_finer_level =
                hierarchy->getPatchLevel(ln+1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->
                                         getRatioToCoarserLevel());
        }
        
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            
            const double* dx = pgeom->getDx();
            const double cell_vol = dx[0]
#if (NDIM > 1)
                *dx[1]
#if (NDIM > 2)
                *dx[2]
#endif
#endif
                ;
            
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > wgt_data =
                patch->getPatchData(d_wgt_idx);
            
            wgt_data->fillAll(cell_vol);
            
            if (ln < finest_hier_level)
            {
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes();
                     ++i)
                {
                    const SAMRAI::hier::Box<NDIM>& refined_box = refined_region_boxes(i);
                    const SAMRAI::hier::Box<NDIM> intersection = patch_box * refined_box;
                    
                    if (!intersection.empty())
                    {
                        wgt_data->fillAll(0.0, intersection);
                    }
                }
            }
        }
    }
    return;
}// resetHierarchyConfiguration

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the ConvergenceMonitor
/// class.
///
//////////////////////////////////////////////////////////////////////

template class SAMRAI::tbox::Pointer<IBAMR::ConvergenceMonitor>;

///////////////////////////////////////////////////////////////////////////////
