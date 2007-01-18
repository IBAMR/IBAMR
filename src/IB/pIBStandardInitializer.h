#ifndef included_pIBStandardInitializer
#define included_pIBStandardInitializer

// Filename: pIBStandardInitializer.h
// Last modified: <18.Jan.2007 15:38:54 boyce@bigboy.nyconnect.com>
// Created on 18 Jan 2007 by Boyce Griffith (boyce@bigboy.nyconnect.com)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodePosnInitStrategy.h>
#include <ibamr/LagSiloDataWriter.h>
#include <ibamr/Stashable.h>

// C++ STDLIB INCLUDES
#include <map>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class pIBStandardInitialier is a concrete
 * LNodePosnInitStrategy that can intialize the configuration of one
 * or more IB structures from input files.  This initializer class is
 * intended to be used with the penalty immersed boundary (pIB) method
 * for "massive" immersed elastic structures, i.e., for elastic
 * structures whose mass densities are different from that of the
 * background fluid.
 *
 * \todo document input database entries
 */
class pIBStandardInitializer
    : public LNodePosnInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    pIBStandardInitializer(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~pIBStandardInitializer();

    /*!
     * \brief Register a Silo data writer with the pIB initializer
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
    pIBStandardInitializer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    pIBStandardInitializer(
        const pIBStandardInitializer& from);

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
    pIBStandardInitializer& operator=(
        const pIBStandardInitializer& that);

    /*!
     * \brief Read the vertex data from one or more input files.
     */
    void readVertexFiles();

    /*!
     * \brief Read the edge data from one or more input files.
     */
    void readEdgeFiles();

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
     * \return The mass of the specified vertex.
     */
    double getVertexMass(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The anchor tag of the specified vertex.
     */
    bool getVertexAnchorTag(
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
    std::vector<std::vector<std::string> > d_base_filenames;

    /*
     * Vertex information.
     */
    std::vector<std::vector<int> > d_num_vertices, d_vertex_offsets;
    std::vector<std::vector<std::vector<double> > > d_vertex_posns, d_vertex_masses;
    std::vector<std::vector<std::vector<int> > > d_vertex_anchor_tags;

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
    std::vector<std::vector<std::multimap<int,Edge> > > d_edge_map;
    std::vector<std::vector<std::map<Edge,double,EdgeComp> > > d_edge_stiffnesses, d_edge_rest_lengths;

    /*
     * Data required to specify connectivity information for
     * visualization purposes.
     */
    std::vector<int> d_global_index_offset;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "pIBStandardInitializer.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_pIBStandardInitializer
