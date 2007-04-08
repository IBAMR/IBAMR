#ifndef included_IBStandardInitializer
#define included_IBStandardInitializer

// Filename: IBStandardInitializer.h
// Last modified: <07.Apr.2007 18:18:42 griffith@box221.cims.nyu.edu>
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
 * \todo Document input database entries and input file formats.
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
        const bool initial_time,
        LDataManager* const lag_manager);

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
        const bool initial_time,
        LDataManager* const lag_manager);

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
     * \brief Configure the Lagrangian Silo data writer to plot the data
     * associated with the specified level of the locally refined Cartesian
     * grid.
     */
    void initializeLagSiloDataWriter(
        const int level_number);

    /*!
     * \brief Read the vertex data from one or more input files.
     */
    void readVertexFiles();

    /*!
     * \brief Read the spring data from one or more input files.
     */
    void readSpringFiles();

    /*!
     * \brief Read the beam data from one or more input files.
     */
    void readBeamFiles();

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
     * non-null.
     */
    void getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.
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
    SAMRAI::tbox::Pointer<LagSiloDataWriter> d_silo_writer;

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
     * Spring information.
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
    std::vector<std::vector<bool> > d_enable_springs;
    std::vector<std::vector<std::multimap<int,Edge> > > d_spring_edge_map;
    std::vector<std::vector<std::map<Edge,double,EdgeComp> > > d_spring_stiffness, d_spring_rest_length;
    std::vector<std::vector<std::map<Edge,int,EdgeComp> > > d_spring_force_fcn_idx;

    std::vector<std::vector<bool> > d_using_uniform_spring_stiffness;
    std::vector<std::vector<double> > d_uniform_spring_stiffness;

    std::vector<std::vector<bool> > d_using_uniform_spring_rest_length;
    std::vector<std::vector<double> > d_uniform_spring_rest_length;

    std::vector<std::vector<bool> > d_using_uniform_spring_force_fcn_idx;
    std::vector<std::vector<int> > d_uniform_spring_force_fcn_idx;

    /*
     * Beam information.
     */
    typedef std::pair<int,int> Neighbors;
    std::vector<std::vector<bool> > d_enable_beams;
    std::vector<std::vector<std::multimap<int,std::pair<Neighbors,double> > > > d_beam_specs;

    std::vector<std::vector<bool> > d_using_uniform_beam_bend_rigidity;
    std::vector<std::vector<double> > d_uniform_beam_bend_rigidity;

    /*
     * Target point information.
     */
    std::vector<std::vector<bool> > d_enable_target_points;
    std::vector<std::vector<std::vector<double> > > d_target_stiffness;

    std::vector<std::vector<bool> > d_using_uniform_target_stiffness;
    std::vector<std::vector<double> > d_uniform_target_stiffness;

    /*
     * Mass information for the pIB method.
     */
    std::vector<std::vector<bool> > d_enable_bdry_mass;
    std::vector<std::vector<std::vector<double> > > d_bdry_mass, d_bdry_mass_stiffness;

    std::vector<std::vector<bool> > d_using_uniform_bdry_mass;
    std::vector<std::vector<double> > d_uniform_bdry_mass;

    std::vector<std::vector<bool> > d_using_uniform_bdry_mass_stiffness;
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
