#ifndef included_IBStandardInitializer
#define included_IBStandardInitializer

// Filename: IBStandardInitializer.h
// Last modified: <12.Jun.2007 18:54:57 griffith@box221.cims.nyu.edu>
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
 * \brief Class IBStandardInitializer is a concrete LNodeInitStrategy that
 * initializes the configuration of one or more Lagrangian structures from input
 * files.
 *
 * \todo Document input database entries.
 *
 * \note "C-style" indices are used for all input files.
 *
 * <HR>
 *
 * <B>Vertex file format</B>
 *
 * Vertex input files end with the extension <TT>".vertex"</TT> and have the
 * following format for two spatial dimensions: \verbatim

 N                   # number of vertices in the file
 x_0       y_0       # (x,y)-coordinates of vertex 0
 x_1       y_1       # (x,y)-coordinates of vertex 1
 ...
 x_{N-1}   y_{N-1}   # (x,y)-coordinates of vertex N-1
 \endverbatim
 *
 * Vertex input files end with the extension <TT>".vertex"</TT> and have the
 * following format for three spatial dimensions: \verbatim

 N                             # number of vertices in the file
 x_0       y_0       z_0       # (x,y,z)-coordinates of vertex 0
 x_1       y_1       z_1       # (x,y,z)-coordinates of vertex 1
 ...
 x_{N-1}   y_{N-1}   z_{N-1}   # (x,y,z)-coordinates of vertex N-1
 \endverbatim
 *
 * <HR>
 *
 * <B>Spring file format</B>
 *
 * Spring input files end with the extension <TT>".spring"</TT> and have the
 * following format: \verbatim

 M                                            # number of links in the file
 i_0   j_0   kappa_0   length_0   fcn_idx_0   # first vertex index, second vertex index, spring constant, rest length, spring function index
 i_1   j_1   kappa_1   length_1   fcn_idx_1
 i_2   j_2   kappa_2   length_2   fcn_idx_2
 ...
 \endverbatim
 *
 * \note There is no restriction on the number of springs that may be associated
 * with any particular node of the Lagrangian mesh.
 *
 * \note The rest length and force function indices are \em optional values.  If
 * they are not provided, by default the rest length will be set to the value \a
 * 0.0 and the force function index will be set to \a 0.  This corresponds to a
 * linear spring with zero rest length.
 *
 * \note Spring specifications are used by class LagSiloDataWriter to construct
 * unstructured mesh representations of the Lagrangian structures.
 * Consequently, even if your structure does not have any springs, it may be
 * worthwhile to generate a spring input file with all spring constants set to
 * \a 0.0.
 *
 * \note \a min(i,j) is always used as the "master" node index when constructing
 * the corresponding IBSpringForceSpec object.
 *
 * \see IBSpringForceGen
 * \see IBSpringForceSpec
 *
 * <HR>
 *
 * <B> Beam file format</B>
 *
 * Beam input files end with the extension <TT>".beam"</TT> and have the
 * following format: \verbatim

 M                           # number of beams in the file
 i_0   j_0   k_0   kappa_0   # first vertex index, second vertex index, third vertex index, bending rigidity
 i_1   j_1   k_1   kappa_1   # first vertex index, second vertex index, third vertex index, bending rigidity
 i_2   j_2   k_2   kappa_2   # first vertex index, second vertex index, third vertex index, bending rigidity
 ...
 \endverbatim
 *
 * \note There is no restriction on the number of beams that may be associated
 * with any particular node of the Lagrangian mesh.
 *
 * \note For each bending-resistant triple \a(i,j,k), it is neccessary that
 * vertex \a j correspond to an "interior" node, i.e., a node that is not the
 * first or last node in the beam.
 *
 * \note The second vertex index is always used as the "master" node index when
 * constructing the corresponding IBBeamForceSpec object.
 *
 * \see IBBeamForceGen
 * \see IBBeamForceSpec
 *
 * <HR>
 *
 * <B>Target point file format</B>
 *
 * Target point input files end with the extension <TT>".target"</TT> and have
 * the following format: \verbatim

 M               # number of target points in the file
 i_0   kappa_0   # vertex index, penalty spring constant
 i_1   kappa_1
 i_2   kappa_2
 ...
 \endverbatim
 *
 * \note Target points are anchored to their \em initial positions by linear
 * springs with the specified spring constants and with zero resting lengths.
 * Consequently, target points approximately enforce internal Dirichlet boundary
 * conditions.  The penalty parameter provides control over the energetic
 * penalty imposed when the position of the Lagrangian immersed boundary point
 * deviates from that of its specified fixed location.
 *
 * \see IBTargetPointForceGen
 * \see IBTargetPointForceSpec
 *
 * <HR>
 *
 * <B>Mass point file format</B>
 *
 * Mass point input files end with the extension <TT>".mass"</TT> and have the
 * following format: \verbatim

 M                           # number of mass points in the file
 i_0   mass_0   kappa_0      # vertex index, point mass, penalty spring constant
 i_1   mass_1   kappa_1
 i_2   mass_2   kappa_2
 ...
 \endverbatim
 * \note Mass points are anchored to "ghost" massive particles by linear springs
 * with the specified spring constants and with zero resting lengths.  The
 * massive particles are "isolated" and simply move according to Newton's laws.
 * The penalty parameter provides control over the energetic penalty imposed
 * when the position of the Lagrangian immersed boundary point deviates from
 * that of its massive copy.
 *
 * <HR>
 *
 * <B>Instrumentation file format</B>
 *
 * Instrumentation input files (specifying the nodes employed to determine the
 * time-dependent positions of flow meters and pressure gauges) end with the
 * extension <TT>".inst"</TT> and have the following format: \verbatim

 M                                      # number of instrumentation points in the file
 i_0   meter_idx_0   meter_node_idx_0   # vertex index, meter index, node index within meter
 i_1   meter_idx_1   meter_node_idx_1
 i_2   meter_idx_2   meter_node_idx_2
 ...
 \endverbatim
 * \note Flow meters and pressure gauges are constructed out of "rings" of
 * immersed boundary points.  The flow is computed by computing the total
 * velocity flux through a web spanning the perimeter of the flow meter.  The
 * pressure is measured at the centroid of each flow meter.
 *
 * Note that each meter may have a different number of nodes specifying its
 * perimeter; however, the values of meter_node_idx associated with a particular
 * meter must be a continuous range of integers, starting with index 0.  E.g.,
 * the following is a valid input file: \verbatim

 6           # number of instrumentation points in the file
 0   0   0   # perimeter of meter 0 consists of vertices 0, 1, and 2
 1   0   1
 2   0   2
 9   1   0   # perimeter of meter 1 consists of vertices 9, 10, and 11
 10  1   1
 11  1   2
 \endverbatim
 *
 * \see IBInstrumentPanel
 * \see IBInstrumentationSpec
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
    virtual
    ~IBStandardInitializer();

    /*!
     * \brief Register a Silo data writer with the IB initializer object.
     */
    void
    registerLagSiloDataWriter(
        SAMRAI::tbox::Pointer<LagSiloDataWriter> silo_writer);

    /*!
     * \brief Determine whether there are any Lagrangian nodes on the specified
     * patch level.
     *
     * \return A boolean value indicating whether Lagrangian data is associated
     * with the given level in the patch hierarchy.
     */
    virtual bool
    getLevelHasLagrangianData(
        const int level_number,
        const bool can_be_refined) const;

    /*!
     * \brief Determine the number of local nodes on the specified patch level.
     *
     * \return The number of local nodes on the specified level.
     */
    virtual int
    getLocalNodeCountOnPatchLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time);

    /*!
     * \brief Initialize the LNodeIndex and LNodeLevel data needed to specify
     * the configuration of the curvilinear mesh on the patch level.
     *
     * \return The number of local nodes initialized on the patch level.
     */
    virtual int
    initializeDataOnPatchLevel(
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
     * \brief Initialize the LNodeLevel data needed to specify the mass and
     * spring constant data required by the penalty IB method.
     *
     * \return The number of local nodes initialized on the patch level.
     */
    virtual int
    initializeMassDataOnPatchLevel(
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
     * When the patch hierarchy is being constructed at the initial simulation
     * time, it is necessary to instruct the gridding algorithm where to place
     * local refinement in order to accomodate portions of the curvilinear mesh
     * that will reside in any yet-to-be-constructed level(s) of the patch
     * hierarchy.
     */
    virtual void
    tagCellsForInitialRefinement(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index);

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBStandardInitializer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBStandardInitializer(
        const IBStandardInitializer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBStandardInitializer&
    operator=(
        const IBStandardInitializer& that);

    /*!
     * \brief Configure the Lagrangian Silo data writer to plot the data
     * associated with the specified level of the locally refined Cartesian
     * grid.
     */
    void
    initializeLagSiloDataWriter(
        const int level_number);

    /*!
     * \brief Read the vertex data from one or more input files.
     */
    void
    readVertexFiles();

    void
    readVertexFile_h5(
        const int ln,
        const int j);

    void
    readVertexFile_ascii(
        const int ln,
        const int j);

    /*!
     * \brief Read the spring data from one or more input files.
     */
    void
    readSpringFiles();

    void
    readSpringFile_h5(
        const int ln,
        const int j);

    void
    readSpringFile_ascii(
        const int ln,
        const int j);

    /*!
     * \brief Read the beam data from one or more input files.
     */
    void
    readBeamFiles();

    /*!
     * \brief Read the target point data from one or more input files.
     */
    void
    readTargetPointFiles();

    /*!
     * \brief Read the boundary mass data from one or more input files.
     */
    void
    readBoundaryMassFiles();

    /*!
     * \brief Read the instrumentation data from one or more input files.
     */
    void
    readInstrumentationFiles();

    /*!
     * \brief Determine the indices of any vertices initially located within the
     * specified patch.
     */
    void
    getPatchVertices(
        std::vector<std::pair<int,int> >& point_indices,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        const int level_number,
        const bool can_be_refined) const;

    /*!
     * \return The cannonical Lagrangian index of the specified vertex.
     */
    int
    getCannonicalLagrangianIndex(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The initial position of the specified vertex.
     */
    std::vector<double>
    getVertexPosn(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The target point penalty force spring constant associated with a
     * particular node.  (Note that if this value is zero for any particular
     * node, there will be no target point penalty force at that node.)
     */
    double
    getVertexTargetStiffness(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The mass associated with a particular node.
     */
    double
    getVertexMass(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The mass spring constant associated with a particular node.
     */
    double
    getVertexMassStiffness(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The instrumentation indices associated with a particular node (or
     * std::make_pair(-1,-1) if there is no instrumentation data associated with
     * that node).
     */
    std::pair<int,int>
    getVertexInstrumentationIndices(
        const std::pair<int,int>& point_index,
        const int level_number) const;

    /*!
     * \return The specification objects associated with the specified vertex.
     */
    std::vector<SAMRAI::tbox::Pointer<Stashable> >
    initializeSpecs(
        const std::pair<int,int>& point_index,
        const int global_index_offset,
        const int level_number) const;

    /*!
     * Read input values, indicated above, from given database.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

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
    SAMRAI::tbox::Pointer<LagSiloDataWriter> d_silo_writer;

    /*
     * The base filenames of the structures are used to generate unique names
     * when registering data with the Silo data writer.
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
        inline bool
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
     * Instrumentation information.
     */
    std::vector<std::vector<bool> > d_enable_instrumentation;
    std::vector<std::vector<std::map<int,std::pair<int,int> > > > d_instrument_idx;

    /*
     * Data required to specify connectivity information for visualization
     * purposes.
     */
    std::vector<int> d_global_index_offset;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "IBStandardInitializer.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBStandardInitializer
