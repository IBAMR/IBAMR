#ifndef included_IBStandardInitializer
#define included_IBStandardInitializer

// Filename: IBStandardInitializer.h
// Last modified: <23.Jan.2007 02:55:30 boyce@bigboy.nyconnect.com>
// Created on 22 Nov 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeInitStrategy.h>
#include <ibamr/LagSiloDataWriter.h>
#include <ibamr/Stashable.h>

// C++ STDLIB INCLUDES
#include <map>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBStandardInitialier is a concrete LNodeInitStrategy
 * that can intialize the initial configuration of one or more IB
 * structures from input files.
 *
 * \todo document input database entries
 */
class IBStandardInitializer
    : public LNodeInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBStandardInitializer(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~IBStandardInitializer();

    /*!
     * \brief Register a Silo data writer with the IB initializer
     * object.
     */
    void registerLagSiloDataWriter(
        SAMRAI::tbox::Pointer<LagSiloDataWriter> silo_writer);

    /*!
     * \brief Determine whether there are any Lagrangian nodes on the
     * specified patch level.
     *
     * \return A boolean value indicating whether Lagrangian data is
     * associated with the given level in the patch hierarchy.
     */
    virtual bool getLevelHasLagrangianData(
        const int level_number,
        const bool can_be_refined) const;

    /*!
     * \brief Determine the number of local nodes on the specified
     * patch level.
     *
     * \return The number of local nodes on the specified level.
     */
    virtual int getLocalNodeCountOnPatchLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * \brief Initialize the LNodeIndex and LNodeLevel data needed to
     * specify the configuration of the curvilinear mesh on the patch
     * level.
     *
     * \return The number of local nodes initialized on the patch
     * level.
     */
    virtual int initializeDataOnPatchLevel(
        const int lag_node_index_idx,
        const int global_index_offset,
        const int local_index_offset,
        SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * @brief Initialize the LNodeLevel data needed to specify the
     * mass and spring constant data required by the penalty IB
     * method.
     *
     * \return The number of local nodes initialized on the patch
     * level.
     */
    virtual int initializeMassDataOnPatchLevel(
        const int global_index_offset,
        const int local_index_offset,
        SAMRAI::tbox::Pointer<LNodeLevelData>& M_data,
        SAMRAI::tbox::Pointer<LNodeLevelData>& K_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * \brief Tag cells for initial refinement.
     *
     * When the patch hierarchy is being constructed at the initial
     * simulation time, it is necessary to instruct the gridding
     * algorithm where to place local refinement in order to
     * accomodate portions of the curvilinear mesh that will reside in
     * any yet-to-be-constructed level(s) of the patch hierarchy.
     */
    virtual void tagCellsForInitialRefinement(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    IBStandardInitializer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    IBStandardInitializer(
        const IBStandardInitializer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBStandardInitializer& operator=(
        const IBStandardInitializer& that);

    /*!
     * \brief Read the vertex data from one or more input files.
     */
    void readVertexFiles();

    /*!
     * \brief Read the edge data from one or more input files.
     */
    void readEdgeFiles();

    /*!
     * \brief Read the target point data from one or more input files.
     */
    void readTargetPointFiles();

    /*!
     * \brief Read the boundary mass data from one or more input
     * files.
     */
    void readBoundaryMassFiles();

    /*!
     * \brief Determine the indices of any vertices initially located
     * within the specified patch.
     */
    void getPatchVertices(
        std::vector<std::pair<int,int> >& point_indices,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        const int level_number,
        const bool can_be_refined) const;

    /*!
     * \return The cannonical Lagrangian index of the specified
     * vertex.
     */
    int getCannonicalLagrangianIndex(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The initial position of the specified vertex.
     */
    std::vector<double> getVertexPosn(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The target point penalty force spring constant
     * associated with a particular node.  (Note that if this value is
     * zero for any particular node, there will be no target point
     * penalty force at that node.)
     */
    double getVertexTargetStiffness(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The mass associated with a particular node.
     */
    double getVertexMass(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The mass spring constant associated with a particular
     * node.
     */
    double getVertexMassStiffness(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The force specification objects associated with the
     * specified vertex.
     */
    std::vector<SAMRAI::tbox::Pointer<Stashable> > initializeForceSpec(
        const std::pair<int,int>& point_index,
        const int global_index_offset,
        const int level_number) const;

    /*!
     * Read input values, indicated above, from given database.
     *
     * When assertion checking is active, the database pointer must be
     * non-null.  Otherwise, all your base are belong to us.
     */
    void getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * The maximum number of levels in the Cartesian grid patch
     * hierarchy.
     */
    int d_max_levels;

    /*
     * The base filenames of the structures are used to generate
     * unique names when registering data with the Silo data writer.
     */
    std::vector<std::vector<std::string> > d_base_filename;

    /*
     * Vertex information.
     */
    std::vector<std::vector<int> > d_num_vertex, d_vertex_offset;
    std::vector<std::vector<std::vector<double> > > d_vertex_posn;

    /*
     * Edge information.
     */
    typedef std::pair<int,int> Edge;
    struct EdgeComp
        : public std::binary_function<Edge,Edge,bool>
    {
        bool
        operator()(
            const Edge& e1,
            const Edge& e2) const
            {
                return (e1.first < e2.first) || (e1.first == e2.first && e1.second < e2.second);
            }
    };
    std::vector<std::vector<bool> > d_enable_edges;
    std::vector<std::vector<std::multimap<int,Edge> > > d_edge_map;
    std::vector<std::vector<std::map<Edge,double,EdgeComp> > > d_edge_stiffness, d_edge_rest_length;

    std::vector<std::vector<bool> > d_use_uniform_edge_stiffness;
    std::vector<std::vector<double> > d_uniform_edge_stiffness;

    std::vector<std::vector<bool> > d_use_uniform_edge_rest_length;
    std::vector<std::vector<double> > d_uniform_edge_rest_length;

    /*
     * Target point information.
     */
    std::vector<std::vector<bool> > d_enable_target_points;
    std::vector<std::vector<std::vector<double> > > d_target_stiffness;

    std::vector<std::vector<bool> > d_use_uniform_target_stiffness;
    std::vector<std::vector<double> > d_uniform_target_stiffness;

    /*
     * Mass information for the pIB method.
     */
    std::vector<std::vector<bool> > d_enable_bdry_mass;
    std::vector<std::vector<std::vector<double> > > d_bdry_mass, d_bdry_mass_stiffness;

    std::vector<std::vector<bool> > d_use_uniform_bdry_mass;
    std::vector<std::vector<double> > d_uniform_bdry_mass;

    std::vector<std::vector<bool> > d_use_uniform_bdry_mass_stiffness;
    std::vector<std::vector<double> > d_uniform_bdry_mass_stiffness;

    /*
     * Data required to specify connectivity information for
     * visualization purposes.
     */
    std::vector<int> d_global_index_offset;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "IBStandardInitializer.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBStandardInitializer
