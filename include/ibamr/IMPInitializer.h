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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_IMPInitializer
#define included_IBAMR_IMPInitializer

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include "ibtk/LInitStrategy.h"
#include "ibtk/LSiloDataWriter.h"

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"

#include "libmesh/id_types.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
class LData;
class LDataManager;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
namespace libMesh
{
class MeshBase;
class Point;
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IMPInitializer is a concrete LInitStrategy that initializes the
 * configuration of one or more Lagrangian structures that are described using
 * the immersed material point method from FE meshes.
 */
class IMPInitializer : public IBTK::LInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IMPInitializer(std::string object_name,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                   SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * \brief Destructor.
     */
    ~IMPInitializer();

    /*!
     * \brief Register a Mesh object with the IB initializer object.
     */
    void registerMesh(libMesh::MeshBase* mesh, int level_number = IBTK::invalid_level_number);

    /*!
     * \brief Register a Silo data writer with the IB initializer object.
     */
    void registerLSiloDataWriter(SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> silo_writer);

    /*!
     * \brief Determine whether there are any Lagrangian nodes on the specified
     * patch level.
     *
     * \return A boolean value indicating whether Lagrangian data is associated
     * with the given level in the patch hierarchy.
     */
    bool getLevelHasLagrangianData(int level_number, bool can_be_refined) const override;

    /*!
     * \return A boolean value indicating whether or not all Lagrangian data is
     * within the computational domain specified by the patch hierarchy.
     */
    bool
    getIsAllLagrangianDataInDomain(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy) const override;

    /*!
     * \return Determine the number of global nodes on the patch level.
     *
     * \return The number of global nodes on the specified level.
     */
    unsigned int
    computeGlobalNodeCountOnPatchLevel(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       int level_number,
                                       double init_data_time,
                                       bool can_be_refined,
                                       bool initial_time) override;

    /*!
     * \brief Determine the number of local nodes on the specified patch level.
     *
     * \return The number of local nodes on the specified level.
     */
    unsigned int computeLocalNodeCountOnPatchLevel(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                                   int level_number,
                                                   double init_data_time,
                                                   bool can_be_refined,
                                                   bool initial_time) override;

    /*!
     * \brief Initialize the structure indexing information on the patch level.
     */
    void initializeStructureIndexingOnPatchLevel(std::map<int, std::string>& strct_id_to_strct_name_map,
                                                 std::map<int, std::pair<int, int> >& strct_id_to_lag_idx_range_map,
                                                 int level_number,
                                                 double init_data_time,
                                                 bool can_be_refined,
                                                 bool initial_time,
                                                 IBTK::LDataManager* l_data_manager) override;

    /*!
     * \brief Initialize the LNode and LData data needed to specify the
     * configuration of the curvilinear mesh on the patch level.
     *
     * \return The number of local nodes initialized on the patch level.
     */
    unsigned int initializeDataOnPatchLevel(int lag_node_index_idx,
                                            unsigned int global_index_offset,
                                            unsigned int local_index_offset,
                                            SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                                            SAMRAI::tbox::Pointer<IBTK::LData> U_data,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                            int level_number,
                                            double init_data_time,
                                            bool can_be_refined,
                                            bool initial_time,
                                            IBTK::LDataManager* l_data_manager) override;

    /*!
     * \brief Tag cells for initial refinement.
     *
     * When the patch hierarchy is being constructed at the initial simulation
     * time, it is necessary to instruct the gridding algorithm where to place
     * local refinement in order to accommodate portions of the curvilinear mesh
     * that will reside in any yet-to-be-constructed level(s) of the patch
     * hierarchy.
     */
    void tagCellsForInitialRefinement(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
                                      double error_data_time,
                                      int tag_index) override;

    /*!
     * \brief Write vertex position in a file.
     *
     * \param filename Name of the file to which position of the vertices is
     * to be written.
     *
     * \param mesh_no Integral order in which various meshes are registered
     * with the class on level \em level_number. Indexing starts from 0.
     *
     * \param level_number Level on which the mesh resides. -1 indicates finest
     * grid level.
     */
    void writeVertexFile(std::string filename, int mesh_no, int level_number = IBTK::invalid_level_number);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IMPInitializer() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IMPInitializer(const IMPInitializer& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IMPInitializer& operator=(const IMPInitializer& that) = delete;

    /*!
     * \brief Configure the Lagrangian Silo data writer to plot the data
     * associated with the specified level of the locally refined Cartesian
     * grid.
     */
    void initializeLSiloDataWriter(int level_number);

    /*!
     * \brief Determine the indices of any vertices initially located within the
     * specified patch.
     */
    void getPatchVertices(std::vector<std::pair<int, int> >& point_indices,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                          int level_number,
                          bool can_be_refined) const;

    /*!
     * \return The canonical Lagrangian index of the specified vertex.
     */
    int getCanonicalLagrangianIndex(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The initial position of the specified vertex.
     */
    const libMesh::Point& getVertexPosn(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * Read input values, indicated above, from given database.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * Patch hierarchy on which we are setting up data and corresponding
     * gridding algorithm.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;
    std::vector<bool> d_level_is_initialized;

    /*
     * Assignment of meshes to level numbers.
     */
    std::vector<std::vector<libMesh::MeshBase*> > d_meshes;

    /*
     * Material point data.
     */
    std::vector<std::vector<int> > d_num_vertex, d_vertex_offset;
    std::vector<std::vector<std::vector<libMesh::Point> > > d_vertex_posn;
    std::vector<std::vector<std::vector<double> > > d_vertex_wgt;
    std::vector<std::vector<std::vector<libMesh::subdomain_id_type> > > d_vertex_subdomain_id;

    /*
     * An (optional) Lagrangian Silo data writer.
     */
    SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> d_silo_writer;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBAMR_HAVE_LIBMESH
#endif //#ifndef included_IBAMR_IMPInitializer
