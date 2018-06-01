// Filename: IBRedundantInitializer.h
// Created on 24 May 2018 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBAMR_IBRedundantInitializer
#define included_IBAMR_IBRedundantInitializer

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "IntVector.h"
#include "boost/array.hpp"
#include "ibamr/IBRodForceSpec.h"
#include "ibtk/LInitStrategy.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Pointer.h"
#include <boost/concept_check.hpp>

namespace IBTK
{
class LData;
class LDataManager;
class Streamable;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBRedundantInitializer is an abstract LInitStrategy that
 * initializes the configuration of one or more Lagrangian structures from input
 * files.
 *
 * \todo Document input database entries.
 *
 */
class IBRedundantInitializer : public IBTK::LInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBRedundantInitializer(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~IBRedundantInitializer();

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
    bool getLevelHasLagrangianData(int level_number, bool can_be_refined) const;

    /*!
     * \brief Determine the number of global nodes on the specified patch level.
     *
     * \return The number of global nodes on the specified level.
     */
    unsigned int
    computeGlobalNodeCountOnPatchLevel(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       int level_number,
                                       double init_data_time,
                                       bool can_be_refined,
                                       bool initial_time);

    /*!
     * \brief Determine the number of local nodes on the specified patch level.
     *
     * \return The number of local nodes on the specified level.
     */
    unsigned int computeLocalNodeCountOnPatchLevel(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                                   int level_number,
                                                   double init_data_time,
                                                   bool can_be_refined,
                                                   bool initial_time);

    /*!
     * \brief Initialize structure specific configurations.
     */
    virtual void init();

    /*!
     * \brief Initialize vertex data programmatically.
     */
    void initializeStructurePosition();

    typedef void (*InitStructureOnLevel)(const unsigned int& strct_num,
                                         const int& level_num,
                                         int& num_vertices,
                                         std::vector<IBTK::Point>& vertex_posn);

    void registerInitStructureFunction(InitStructureOnLevel fcn);

    /*!
     * \brief Initialize spring data programmatically.
     */
    void initializeSprings();

    /*
     * Edge data structures.
     */
    typedef std::pair<int, int> Edge;
    struct EdgeComp : public std::binary_function<Edge, Edge, bool>
    {
        inline bool operator()(const Edge& e1, const Edge& e2) const
        {
            return (e1.first < e2.first) || (e1.first == e2.first && e1.second < e2.second);
        }
    };

    struct SpringSpec
    {
        std::vector<double> parameters;
        int force_fcn_idx;
    };
    typedef void (*InitSpringDataOnLevel)(const unsigned int& strct_num,
                                          const int& level_num,
                                          std::multimap<int, Edge>& spring_map,
                                          std::map<Edge, SpringSpec, EdgeComp>& spring_spec);

    void registerInitSpringDataFunction(InitSpringDataOnLevel fcn);

    /*
     * \brief Initialize xspring data programmatically.
     */
    void initializeXSprings();

    struct XSpringSpec
    {
        std::vector<double> parameters;
        int force_fcn_idx;
    };

    typedef void (*InitXSpringDataOnLevel)(const unsigned int& strct_num,
                                           const int& level_num,
                                           std::multimap<int, Edge>& xspring_map,
                                           std::map<Edge, XSpringSpec, EdgeComp> xspring_spec);

    void registerInitXSpringDataFunction(InitXSpringDataOnLevel fcn);

    /*!
     * \brief Initialize beam data programmatically.
     */
    void initializeBeams();

    struct BeamSpec
    {
        std::pair<int, int> neighbor_idxs;
        double bend_rigidity;
        IBTK::Vector curvature;
    };
    typedef void (*InitBeamDataOnLevel)(const unsigned int& strct_num,
                                        const int& level_num,
                                        std::multimap<int, BeamSpec>& beam_spec);

    void registerInitBeamDataFunction(InitBeamDataOnLevel fcn);

    /*!
     * \brief Initialize director and rod data programmatically.
     */
    void initializeDirectorAndRods();

    struct RodSpec
    {
        boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> properties;
    };
    typedef void (*InitDirectorAndRodOnLevel)(const unsigned int& strct_num,
                                              const int& level_num,
                                              std::vector<std::vector<double> >& director_spec,
                                              std::multimap<int, Edge>& rod_edge_map,
                                              std::map<Edge, RodSpec, EdgeComp>& rod_spec);

    void registerInitDirectorAndRodFunction(InitDirectorAndRodOnLevel fcn);

    /*!
     * \brief Initialize massive point data programmatically.
     */
    void initializeBoundaryMass();

    struct BdryMassSpec
    {
        double bdry_mass, stiffness;
    };

    typedef void (*InitBoundaryMassOnLevel)(const unsigned int& strct_num,
                                            const int& level_num,
                                            std::vector<BdryMassSpec>& bdry_mass_spec);

    void registerInitBoundaryMassFunction(InitBoundaryMassOnLevel fcn);

    /*!
     * \brief Initialize target point data programmatically.
     */
    void initializeTargetPts();

    struct TargetSpec
    {
        double stiffness, damping;
    };

    typedef void (*InitTargetPtOnLevel)(const unsigned int& strct_num,
                                        const int& level_num,
                                        std::vector<TargetSpec>& tg_pt_spec);

    void registerInitTargetPtFunction(InitTargetPtOnLevel fcn);

    struct AnchorSpec
    {
        bool is_anchor_point;
    };

    typedef void (*InitAnchorPtOnLevel)(const unsigned int& strct_num,
                                        const int& level_num,
                                        std::vector<AnchorSpec>& anchor_pt_spec);

    void registerInitAnchorPtFunction(InitAnchorPtOnLevel fcn);

    void initializeAnchorPts();

    /*!
     * \brief Initialize instrumentation data.
     *
     * \note Instruments can not currently be implemented using IBRedundantInitializer
     */
    void initializeInstrumentationData();

    /*!
     * \brief Initialize source/sink data.
     *
     * \note Sources/sinks can not currently be implemented using IBRedundantInitializer
     */
    void initializeSourceData();

    /*!
     * \brief Initialize the structure indexing information on the patch level.
     */
    void initializeStructureIndexingOnPatchLevel(std::map<int, std::string>& strct_id_to_strct_name_map,
                                                 std::map<int, std::pair<int, int> >& strct_id_to_lag_idx_range_map,
                                                 int level_number,
                                                 double init_data_time,
                                                 bool can_be_refined,
                                                 bool initial_time,
                                                 IBTK::LDataManager* l_data_manager);

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
                                            IBTK::LDataManager* l_data_manager);

    /*!
     * \brief Initialize the LData needed to specify the mass and spring
     * constant data required by the penalty IB method.
     *
     * \return The number of local nodes initialized on the patch level.
     */
    unsigned int initializeMassDataOnPatchLevel(unsigned int global_index_offset,
                                                unsigned int local_index_offset,
                                                SAMRAI::tbox::Pointer<IBTK::LData> M_data,
                                                SAMRAI::tbox::Pointer<IBTK::LData> K_data,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                                int level_number,
                                                double init_data_time,
                                                bool can_be_refined,
                                                bool initial_time,
                                                IBTK::LDataManager* l_data_manager);

    /*!
     * \brief Initialize the LNode data needed to specify director vectors
     * required by some material models.
     *
     * \return The number of local nodes initialized on the patch level.
     */
    unsigned int
    initializeDirectorDataOnPatchLevel(unsigned int global_index_offset,
                                       unsigned int local_index_offset,
                                       SAMRAI::tbox::Pointer<IBTK::LData> D_data,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       int level_number,
                                       double init_data_time,
                                       bool can_be_refined,
                                       bool initial_time,
                                       IBTK::LDataManager* l_data_manager);

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
                                      int tag_index);

    /*!
     * \brief Set the names of the structures on a given level.
     *
     * \note The structure will be initialized in the same order as the supplied vector.
     */
    void setStructureNamesOnLevel(const int& level_num, const std::vector<std::string>& strct_names);

protected:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBRedundantInitializer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBRedundantInitializer(const IBRedundantInitializer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBRedundantInitializer& operator=(const IBRedundantInitializer& that);

    /*!
     * \brief Configure the Lagrangian Silo data writer to plot the data
     * associated with the specified level of the locally refined Cartesian
     * grid.
     */
    void initializeLSiloDataWriter(int level_number);

    /*!
     * \brief Determine the indices of any vertices initially owned by the
     * specified patch.
     */
    void getPatchVertices(std::vector<std::pair<int, int> >& point_indices,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy) const;

    /*!
     * \brief Determine the indices of any vertices associated with a given
     * level number initially located within the specified patch.
     */
    void getPatchVerticesAtLevel(std::vector<std::pair<int, int> >& point_indices,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 int level_number) const;

    /*!
     * \return The canonical Lagrangian index of the specified vertex.
     */
    int getCanonicalLagrangianIndex(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The initial position of the specified vertex.
     */
    IBTK::Point getVertexPosn(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The initial position of the specified vertex.
     */
    IBTK::Point getShiftedVertexPosn(const std::pair<int, int>& point_index,
                                     int level_number,
                                     const double* domain_x_lower,
                                     const double* domain_x_upper,
                                     const SAMRAI::hier::IntVector<NDIM>& periodic_shift) const;

    /*!
     * \return The target point specifications associated with a particular
     * node.
     */
    const TargetSpec& getVertexTargetSpec(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The anchor point specifications associated with a particular
     * node.
     */
    const AnchorSpec& getVertexAnchorSpec(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The massive boundary point specifications associated with a
     * particular node.
     */
    const BdryMassSpec& getVertexBdryMassSpec(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The directors associated with a particular node.
     */
    const std::vector<double>& getVertexDirectors(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The instrumentation indices associated with a particular node (or
     * std::make_pair(-1,-1) if there is no instrumentation data associated with
     * that node).
     */
    std::pair<int, int> getVertexInstrumentationIndices(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The source indices associated with a particular node (or -1 if
     * there is no source data associated with that node).
     */
    int getVertexSourceIndices(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The specification objects associated with the specified vertex.
     */
    virtual std::vector<SAMRAI::tbox::Pointer<IBTK::Streamable> >
    initializeNodeData(const std::pair<int, int>& point_index,
                       unsigned int global_index_offset,
                       int level_number) const;

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
     * The maximum number of levels in the Cartesian grid patch hierarchy and a
     * vector of boolean values indicating whether a particular level has been
     * initialized yet.
     */
    int d_max_levels;
    std::vector<bool> d_level_is_initialized;

    /*
     * An (optional) Lagrangian Silo data writer.
     */
    SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> d_silo_writer;

    /*
     * The base filenames of the structures are used to generate unique names
     * when registering data with the Silo data writer.
     */
    std::vector<std::vector<std::string> > d_base_filename;

    /*
     * Optional shift and scale factors.
     *
     * \note These shift and scale factors are applied to ALL structures read in
     * by this reader.
     *
     * \note The scale factor is applied both to positions and to spring rest
     * lengths.
     *
     * \note The shift factor should have the same units as the positions in the
     * input files, i.e., X_final = scale*(X_initial + shift).
     */
    double d_length_scale_factor;
    IBTK::Vector d_posn_shift;

    /*
     * Vertex information.
     */
    std::vector<std::vector<int> > d_num_vertex, d_vertex_offset;
    std::vector<std::vector<std::vector<IBTK::Point> > > d_vertex_posn;

    /*
     * Spring information.
     */
    std::vector<std::vector<std::multimap<int, Edge> > > d_spring_edge_map;

    std::vector<std::vector<std::map<Edge, SpringSpec, EdgeComp> > > d_spring_spec_data;

    /*
     * Crosslink spring ("x-spring") information.
     */
    std::vector<std::vector<std::multimap<int, Edge> > > d_xspring_edge_map;

    std::vector<std::vector<std::map<Edge, XSpringSpec, EdgeComp> > > d_xspring_spec_data;

    /*
     * Beam information.
     */
    std::vector<std::vector<std::multimap<int, BeamSpec> > > d_beam_spec_data;

    /*
     * Rod information.
     */
    std::vector<std::vector<std::multimap<int, Edge> > > d_rod_edge_map;

    std::vector<std::vector<std::map<Edge, RodSpec, EdgeComp> > > d_rod_spec_data;

    /*
     * Target point information.
     */
    std::vector<std::vector<std::vector<TargetSpec> > > d_target_spec_data;

    /*
     * Anchor point information.
     */
    std::vector<std::vector<std::vector<AnchorSpec> > > d_anchor_spec_data;

    /*
     * Mass information for the pIB method.
     */
    std::vector<std::vector<std::vector<BdryMassSpec> > > d_bdry_mass_spec_data;

    /*
     * Orthonormal directors for the generalized IB method.
     */
    std::vector<std::vector<std::vector<std::vector<double> > > > d_directors;

    /*
     * Instrumentation information.
     */
    std::vector<std::vector<std::map<int, std::pair<int, int> > > > d_instrument_idx;

    /*
     * Source information.
     */
    std::vector<std::vector<std::map<int, int> > > d_source_idx;

    /*
     * Data required to specify connectivity information for visualization
     * purposes.
     */
    std::vector<unsigned int> d_global_index_offset;

    /*!
     * Check if user defined data has been processed.
     */
    bool d_data_processed;

private:
    /*
     * Functions used to initialize structures programmatically.
     */
    InitStructureOnLevel d_init_structure_on_level_fcn;
    InitSpringDataOnLevel d_init_spring_on_level_fcn;
    InitXSpringDataOnLevel d_init_xspring_on_level_fcn;
    InitBeamDataOnLevel d_init_beam_on_level_fcn;
    InitDirectorAndRodOnLevel d_init_director_and_rod_on_level_fcn;
    InitBoundaryMassOnLevel d_init_boundary_mass_on_level_fcn;
    InitTargetPtOnLevel d_init_target_pt_on_level_fcn;
    InitAnchorPtOnLevel d_init_anchor_pt_on_level_fcn;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBRedundantInitializer
