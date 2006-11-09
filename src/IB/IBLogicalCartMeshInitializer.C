// Filename: IBLogicalCartMeshInitializer.C
// Last modified: <09.Nov.2006 00:15:55 boyce@bigboy.nyconnect.com>
// Created on 06 Dec 2005 by Boyce Griffith (boyce@boyce.cims.nyu.edu).

#include "IBLogicalCartMeshInitializer.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/SpringForceSpec.h>
#include <ibamr/LNodeIndexData.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <Index.h>
#include <tbox/MPI.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <fstream>
#include <iostream>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBLogicalCartMeshInitializer::IBLogicalCartMeshInitializer(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name),
      d_mesh_name(),
      d_mesh_file(),
      d_ncomp(),
      d_mesh_offset(),
      d_mesh_rank(),
      d_mesh_level(),
      d_mesh_dim(),
      d_mesh_periodic()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
#endif

    // Register the SpringForceSpec object with the StashableManager
    // class.
    SpringForceSpec::registerWithStashableManager();

    // Get the input filename.
    string input_filename;
    if (input_db->keyExists("input_filename"))
    {
        input_filename = input_db->getString("input_filename");
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Key data `input_filename' not found in input.");
    }

    // XXXX dummy data
    d_mesh_name.insert("dummy");
    d_mesh_file["dummy"] = "dummy";
    d_ncomp["dummy"] = 1;
    d_mesh_rank["dummy"] = std::vector<int>(1,2);
    d_mesh_level["dummy"] = std::vector<int>(1,-1);
    d_mesh_dim["dummy"] = std::vector<std::vector<int> >(1, std::vector<int>(2,10));
    d_mesh_periodic["dummy"] = std::vector<std::vector<bool> >(1, std::vector<bool>(2,false));

    // Ensure that the rank of each component is at most NDIM.
    for (std::map<std::string,std::vector<int> >::const_iterator map_it = d_mesh_rank.begin();
         map_it != d_mesh_rank.end(); ++map_it)
    {
        for (std::vector<int>::const_iterator it = (*map_it).second.begin();
             it != (*map_it).second.end(); ++it)
        {
            if ((*it) > NDIM || (*it) < 1)
            {
                TBOX_ERROR(d_object_name << "::IBLogicalCartMeshInitializer():\n"
                           << "  the rank of each component in a logically Cartesian mesh must be\n"
                           << "  between 1 and " << NDIM << "\n");
            }
        }
    }

    // Compute the "mesh offset" for each component.
    int offset = 0;
    for (std::set<std::string>::const_iterator set_it = d_mesh_name.begin();
         set_it != d_mesh_name.end(); ++set_it)
    {
        const std::string& name = *set_it;
        d_mesh_offset[name] = offset;
        offset += accumulate(d_mesh_rank[name].begin(),
                             d_mesh_rank[name].end(),
                             1, multiplies<int>());
    }

    return;
}// IBLogicalCartMeshInitializer

IBLogicalCartMeshInitializer::~IBLogicalCartMeshInitializer()
{
    return;
}// ~IBLogicalCartMeshInitializer

namespace
{
// When the level requested for a particular component is negative
// or greater than the finest level number in the hierarchy, that
// component is assigned to the finest level of the patch
// hierarchy.
inline bool is_assigned_level(
    const int requested_level,
    const int level_number,
    const bool can_be_refined)
{
    return (( level_number == requested_level ) ||
            ( (requested_level < 0 || requested_level > level_number)
              && !can_be_refined ));
}
}

bool
IBLogicalCartMeshInitializer::getLevelHasLagrangianData(
    const int level_number,
    const bool can_be_refined) const
{
    // Check to see if any of the components of any of the meshes are
    // assigned to the present level.
    //
    for (std::map<std::string,std::vector<int> >::const_iterator map_it = d_mesh_level.begin();
         map_it != d_mesh_level.end(); ++map_it)
    {
        for (std::vector<int>::const_iterator it = (*map_it).second.begin();
             it != (*map_it).second.end(); ++it)
        {
            const int requested_level = (*it);
            if (is_assigned_level(requested_level, level_number, can_be_refined))
            {
                return true;
            }
        }
    }
    return false;
}// getLevelHasLagrangianData

int
IBLogicalCartMeshInitializer::getLocalNodeCountOnPatchLevel(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    // Count the number of local mesh nodes whose initial positions
    // lie in the local portion of the specified level of the patch
    // hierarchy.
    int num_nodes = 0;
    for (std::map<std::string,std::vector<int> >::const_iterator map_it = d_mesh_level.begin();
         map_it != d_mesh_level.end(); ++map_it)
    {
        for (std::vector<int>::const_iterator it = (*map_it).second.begin();
             it != (*map_it).second.end(); ++it)
        {
            const int requested_level = (*it);
            if (is_assigned_level(requested_level, level_number, can_be_refined))
            {
                const std::string& name = (*map_it).first;
                const int comp = it - (*map_it).second.begin();
                num_nodes += getMeshLocalNodeCount(
                    name, comp, hierarchy->getPatchLevel(level_number));
            }
        }
    }
    return num_nodes;
}// getLocalNodeCountOnPatchLevel

void
IBLogicalCartMeshInitializer::initializeDataOnPatchLevel(
    const int lag_node_index_idx,
    SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    // Initialize the local data on the specified level of the patch
    // hierarchy.
    for (std::map<std::string,std::vector<int> >::const_iterator map_it = d_mesh_level.begin();
         map_it != d_mesh_level.end(); ++map_it)
    {
        for (std::vector<int>::const_iterator it = (*map_it).second.begin();
             it != (*map_it).second.end(); ++it)
        {
            const int requested_level = (*it);
            if (is_assigned_level(requested_level, level_number, can_be_refined))
            {
                const std::string& name = (*map_it).first;
                const int comp = it - (*map_it).second.begin();
                initializeLocalData(
                    name, comp, lag_node_index_idx, X_data,
                    hierarchy->getPatchLevel(level_number));
            }
        }
    }
    return;
}// initializeDataOnPatchLevel

void
IBLogicalCartMeshInitializer::tagCellsForInitialRefinement(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    // Tag cells for refinement on the specified level of the patch
    // hierarchy.
    for (std::map<std::string,std::vector<int> >::const_iterator map_it = d_mesh_level.begin();
         map_it != d_mesh_level.end(); ++map_it)
    {
        for (std::vector<int>::const_iterator it = (*map_it).second.begin();
             it != (*map_it).second.end(); ++it)
        {
            const int requested_level = (*it);
            if (requested_level > level_number || requested_level < 0)
            {
                const std::string& name = (*map_it).first;
                const int comp = it - (*map_it).second.begin();
                tagLocalCellsForRefinement(
                    name, comp, tag_index,
                    hierarchy->getPatchLevel(level_number));
            }
        }
    }
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

int
IBLogicalCartMeshInitializer::getMeshLocalNodeCount(
    const std::string& mesh_name,
    const int comp_num,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& patch_level)
{
    return 0;
}// getMeshLocalNodeCount

void
IBLogicalCartMeshInitializer::initializeLocalData(
    const std::string& mesh_name,
    const int comp_num,
    const int lag_node_index_idx,
    SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& patch_level)
{
    return;
}// initializeLocalData

void
IBLogicalCartMeshInitializer::tagLocalCellsForRefinement(
    const std::string& mesh_name,
    const int comp_num,
    const int tag_index,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >& patch_level)
{
    return;
}// tagLocalCellsForRefinement

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBLogicalCartMeshInitializer>;

//////////////////////////////////////////////////////////////////////////////
