// Filename: IBStandardInitializer.h
// Created on 22 Nov 2006 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#ifndef included_IBStandardInitializer
#define included_IBStandardInitializer

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
 * \brief Class IBStandardInitializer is a concrete LInitStrategy that
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
 * following format for two-dimensional models:
 \verbatim
 N           # number of vertices in the file
 x_0   y_0   # (x,y)-coordinates of vertex 0
 x_1   y_1   # (x,y)-coordinates of vertex 1
 x_2   y_2   # (x,y)-coordinates of vertex 2
 ...
 \endverbatim
 *
 * Vertex input files end with the extension <TT>".vertex"</TT> and have the
 * following format for three-dimensional models:
 \verbatim
 N                 # number of vertices in the file
 x_0   y_0   z_0   # (x,y,z)-coordinates of vertex 0
 x_1   y_1   z_1   # (x,y,z)-coordinates of vertex 1
 x_2   y_2   z_2   # (x,y,z)-coordinates of vertex 2
 ...
 \endverbatim
 *
 * <HR>
 *
 * <B>Spring file format</B>
 *
 * Spring input files end with the extension <TT>".spring"</TT> and have the
 * following format:
 \verbatim
 M                                            # number of links in the file
 i_0   j_0   kappa_0   length_0   fcn_idx_0   # first vertex index, second vertex index, spring
 constant, rest length, spring function index
 i_1   j_1   kappa_1   length_1   fcn_idx_1
 i_2   j_2   kappa_2   length_2   fcn_idx_2
 ...
 \endverbatim
 *
 * \note There is no restriction on the number of springs that may be associated
 * with any particular node of the Lagrangian mesh.
 *
 * \note The rest length and force function index are \em optional values.  If
 * they are not provided, by default the rest length will be set to the value \a
 * 0.0 and the force function index will be set to \a 0.  This corresponds to a
 * linear spring with zero rest length.
 *
 * \note Spring specifications are used by class LSiloDataWriter to construct
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
 * <B>Crosslink spring file format</B>
 *
 * Crosslink spring ("x-spring") input files end with the extension
 * <TT>".xspring"</TT> and have the following format:
 \verbatim
 M                                            # number of links in the file
 i_0   j_0   kappa_0   length_0   fcn_idx_0   # first vertex index, second vertex index, spring
 constant, rest length, spring function index
 i_1   j_1   kappa_1   length_1   fcn_idx_1
 i_2   j_2   kappa_2   length_2   fcn_idx_2
 ...
 \endverbatim
 *
 * \note Unlike standard spring files, in which all indices are required to
 * refer to points within a particular structure, x-spring files may connect
 * points from different structures.  Consequently, the node indices in an
 * x-spring file must be \em global indices.  Notice that global indices are
 * determined by the order in which the structures are specified in the input
 * file.  Changes in the order in which structures are specified necessarily
 * change the global indexing scheme.
 *
 * \note Crosslink springs may connect only structures assigned to the \em same
 * level of the locally refined grid.
 *
 * \note There is no restriction on the number of x-springs that may be
 * associated with any particular node of the Lagrangian mesh.
 *
 * \note The rest length and force function index are \em optional values.  If
 * they are not provided, then by default the rest length will be set to the
 * value \a 0.0 and the force function index will be set to \a 0.  This
 * corresponds to a linear spring with zero rest length.
 *
 * \note Crosslink spring specifications are used by class LSiloDataWriter to
 * construct unstructured mesh representations of the Lagrangian structures.
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
 * following format:
 \verbatim
 M                           # number of beams in the file
 i_0   j_0   k_0   kappa_0   # first vertex index, second vertex index, third vertex index,
 bending
 rigidity
 i_1   j_1   k_1   kappa_1
 i_2   j_2   k_2   kappa_2
 ...
 \endverbatim
 *
 * \note There is no restriction on the number of beams that may be associated
 * with any particular node of the Lagrangian mesh.
 *
 * \note For each bending-resistant triple \a(i,j,k), it is necessary that
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
 * <B> Rod file format</B>
 *
 * Rod input files end with the extension <TT>".rod"</TT> and have the following
 * format:
 \verbatim
 M                                                                                          #
 number
 of rods in the file
 i_0   j_0   ds_0   a1_0   a2_0   a3_0   b1_0   b2_0   b3_0   kappa1_0   kappa2_0   tau_0   #
 first
 vertex index, second vertex index, material parameters
 i_1   j_1   ds_1   a1_1   a2_1   a3_1   b1_1   b2_1   b3_1   kappa1_1   kappa2_1   tau_1
 i_2   j_2   ds_2   a1_2   a2_2   a3_2   b1_2   b2_2   b3_2   kappa1_2   kappa2_2   tau_2
 ...
 \endverbatim
 *
 * \note There is no restriction on the number of rods that may be associated
 * with any particular node of the Lagrangian mesh.
 *
 * \note The first vertex index is always used as the "master" node index when
 * constructing the corresponding IBRodForceSpec object.
 *
 * \note The parameters kappa1, kappa2, and tau (the intrinsic curvatures and
 * twist of the rod) are optional.  If not provided in the input file, they are
 * assumed to be zero.
 *
 * \see IBKirchhoffRodForceGen
 * \see IBRodForceSpec
 *
 * <HR>
 *
 * <B>Target point file format</B>
 *
 * Target point input files end with the extension <TT>".target"</TT> and have
 * the following format:
 \verbatim
 M                       # number of target points in the file
 i_0   kappa_0   eta_0   # vertex index, penalty spring constant, penalty damping coefficient
 i_1   kappa_1   eta_1
 i_2   kappa_2   eta_2
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
 * \note Damping coefficients \f$ \eta \f$ are optional and are set to \a 0.0 if
 * not supplied.  Target points are "anchored" in place using Kelvin-Voigt
 * viscoelastic elements.
 *
 * \see IBTargetPointForceGen
 * \see IBTargetPointForceSpec
 *
 * <HR>
 *
 * <B>Anchor point file format</B>
 *
 * Anchor point input files end with the extension <TT>".anchor"</TT> and have
 * the following format:
 \verbatim
 M                           # number of anchor points in the file
 i_0                         # vertex index
 i_1
 i_2
 ...
 \endverbatim
 * \note Anchor points are immersed boundary nodes that are "anchored" in place.
 * Such points neither spread force nor interpolate velocity.
 *
 * <HR>
 *
 * <B>Mass point file format</B>
 *
 * Mass point input files end with the extension <TT>".mass"</TT> and have the
 * following format:
 \verbatim
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
 * extension <TT>".inst"</TT> and have the following format:
 \verbatim
 M                                      # number of instruments in the file
 meter_name_0                           # meter names
 meter_name_1
 meter_name_2
 ...
 N                                      # number of instrumentation points in the file
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
 * the following is a valid input file:
 *
 \verbatim
 2           # number of instruments in the file
 meter0      # meter names
 meter1
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
 *
 * <HR>
 *
 * <B>Source/sink file format</B>
 *
 * Source/sink input files (specifying the nodes employed to determine the
 * time-dependent positions of internal sources and sinks) end with the
 * extension <TT>".source"</TT> and have the following format:
 \verbatim
 M                    # number of sources/sinks in the file
 source_name_0        # source/sink names
 source_name_1
 source_name_2
 ...
 source_radius_0      # source/sink radii
 source_radius_1
 source_radius_2
 ...
 N                    # number of source/sink points in the file
 i_0   source_idx_0   # vertex index, source/sink index
 i_1   source_idx_1
 i_2   source_idx_2
 ...
 \endverbatim
 * \note The position of each internal source/sink is the arithmetic mean of the
 * positions of the nodes that are associated with that source/sink.
 *
 \verbatim
 2         # number of sources/sinks in the file
 source0   # source/sink names
 source1
 1.0       # source/sink radii
 0.5
 6         # number of source/sink points in the file
 0   0     # position of source0 is determined from vertices 0, 1, and 2
 1   0
 2   0
 9   1     # position of source1 is determined from vertices 9, 10, and 11
 10  1
 11  1
 \endverbatim
 *
 * \see IBStandardSourceGenerator
 * \see IBSourceSpec
 *
 * <HR>
 *
 * <B>Director file format</B>
 *
 * Orthonormal director vector input files end with the extension
 * <TT>".director"</TT> and have the following format, independent of spatial
 * dimension:
 \verbatim
 N                         # number of triads in the file
 D0_x_0   D0_y_0   D0_z_0  # coordinates of director D0 associated with vertex 0
 D1_x_0   D1_y_0   D1_z_0  # coordinates of director D1 associated with vertex 0
 D2_x_0   D2_y_0   D2_z_0  # coordinates of director D2 associated with vertex 0
 D0_x_1   D0_y_1   D0_z_1  # coordinates of director D0 associated with vertex 1
 D1_x_1   D1_y_1   D1_z_1  # coordinates of director D1 associated with vertex 1
 D2_x_1   D2_y_1   D2_z_1  # coordinates of director D2 associated with vertex 1
 D0_x_2   D0_y_2   D0_z_2  # coordinates of director D0 associated with vertex 2
 D1_x_2   D1_y_2   D1_z_2  # coordinates of director D1 associated with vertex 2
 D2_x_2   D2_y_2   D2_z_2  # coordinates of director D2 associated with vertex 2
 ...
 \endverbatim
*/
class IBStandardInitializer : public IBTK::LInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBStandardInitializer(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~IBStandardInitializer();

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
    IBStandardInitializer(const IBStandardInitializer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBStandardInitializer& operator=(const IBStandardInitializer& that);

    /*!
     * \brief Configure the Lagrangian Silo data writer to plot the data
     * associated with the specified level of the locally refined Cartesian
     * grid.
     */
    void initializeLSiloDataWriter(int level_number);

    /*!
     * \brief Read the vertex data from one or more input files.
     */
    void readVertexFiles(const std::string& extension);

    /*!
     * \brief Read the spring data from one or more input files.
     */
    void readSpringFiles(const std::string& file_extension, bool input_uses_global_idxs);

    /*!
     * \brief Read the crosslink spring ("x-spring") data from one or more input
     * files.
     */
    void readXSpringFiles(const std::string& file_extension, bool input_uses_global_idxs);

    /*!
     * \brief Read the beam data from one or more input files.
     */
    void readBeamFiles(const std::string& file_extension, bool input_uses_global_idxs);

    /*!
     * \brief Read the rod data from one or more input files.
     */
    void readRodFiles(const std::string& file_extension, bool input_uses_global_idxs);

    /*!
     * \brief Read the target point data from one or more input files.
     */
    void readTargetPointFiles(const std::string& file_extension);

    /*!
     * \brief Read the anchor point data from one or more input files.
     */
    void readAnchorPointFiles(const std::string& file_extension);

    /*!
     * \brief Read the boundary mass data from one or more input files.
     */
    void readBoundaryMassFiles(const std::string& file_extension);

    /*!
     * \brief Read the director data from one or more input files.
     */
    void readDirectorFiles(const std::string& file_extension);

    /*!
     * \brief Read the instrumentation data from one or more input files.
     */
    void readInstrumentationFiles(const std::string& file_extension);

    /*!
     * \brief Read the source/sink data from one or more input files.
     */
    void readSourceFiles(const std::string& file_extension);

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
    struct TargetSpec;

    const TargetSpec& getVertexTargetSpec(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The anchor point specifications associated with a particular
     * node.
     */
    struct AnchorSpec;

    const AnchorSpec& getVertexAnchorSpec(const std::pair<int, int>& point_index, int level_number) const;

    /*!
     * \return The massive boundary point specifications associated with a
     * particular node.
     */
    struct BdryMassSpec;

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
    std::vector<SAMRAI::tbox::Pointer<IBTK::Streamable> > initializeNodeData(const std::pair<int, int>& point_index,
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
     * The boolean value determines whether file read batons are employed to
     * prevent multiple MPI processes from accessing the same input files
     * simultaneously.
     */
    bool d_use_file_batons;

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

    /*
     * Spring information.
     */
    std::vector<std::vector<bool> > d_enable_springs;

    std::vector<std::vector<std::multimap<int, Edge> > > d_spring_edge_map;

    struct SpringSpec
    {
        std::vector<double> parameters;
        int force_fcn_idx;
    };
    std::vector<std::vector<std::map<Edge, SpringSpec, EdgeComp> > > d_spring_spec_data;

    std::vector<std::vector<bool> > d_using_uniform_spring_stiffness;
    std::vector<std::vector<double> > d_uniform_spring_stiffness;

    std::vector<std::vector<bool> > d_using_uniform_spring_rest_length;
    std::vector<std::vector<double> > d_uniform_spring_rest_length;

    std::vector<std::vector<bool> > d_using_uniform_spring_force_fcn_idx;
    std::vector<std::vector<int> > d_uniform_spring_force_fcn_idx;

    /*
     * Crosslink spring ("x-spring") information.
     */
    std::vector<std::vector<bool> > d_enable_xsprings;

    std::vector<std::vector<std::multimap<int, Edge> > > d_xspring_edge_map;

    struct XSpringSpec
    {
        std::vector<double> parameters;
        int force_fcn_idx;
    };
    std::vector<std::vector<std::map<Edge, XSpringSpec, EdgeComp> > > d_xspring_spec_data;

    std::vector<std::vector<bool> > d_using_uniform_xspring_stiffness;
    std::vector<std::vector<double> > d_uniform_xspring_stiffness;

    std::vector<std::vector<bool> > d_using_uniform_xspring_rest_length;
    std::vector<std::vector<double> > d_uniform_xspring_rest_length;

    std::vector<std::vector<bool> > d_using_uniform_xspring_force_fcn_idx;
    std::vector<std::vector<int> > d_uniform_xspring_force_fcn_idx;

    /*
     * Beam information.
     */
    std::vector<std::vector<bool> > d_enable_beams;

    struct BeamSpec
    {
        std::pair<int, int> neighbor_idxs;
        double bend_rigidity;
        IBTK::Vector curvature;
    };
    std::vector<std::vector<std::multimap<int, BeamSpec> > > d_beam_spec_data;

    std::vector<std::vector<bool> > d_using_uniform_beam_bend_rigidity;
    std::vector<std::vector<double> > d_uniform_beam_bend_rigidity;

    std::vector<std::vector<bool> > d_using_uniform_beam_curvature;
    std::vector<std::vector<IBTK::Vector> > d_uniform_beam_curvature;

    /*
     * Rod information.
     */
    std::vector<std::vector<bool> > d_enable_rods;

    std::vector<std::vector<std::multimap<int, Edge> > > d_rod_edge_map;

    struct RodSpec
    {
        boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> properties;
    };
    std::vector<std::vector<std::map<Edge, RodSpec, EdgeComp> > > d_rod_spec_data;

    std::vector<std::vector<bool> > d_using_uniform_rod_properties;
    std::vector<std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> > > d_uniform_rod_properties;

    /*
     * Target point information.
     */
    std::vector<std::vector<bool> > d_enable_target_points;

    struct TargetSpec
    {
        double stiffness, damping;
    };
    std::vector<std::vector<std::vector<TargetSpec> > > d_target_spec_data;

    std::vector<std::vector<bool> > d_using_uniform_target_stiffness;
    std::vector<std::vector<double> > d_uniform_target_stiffness;

    std::vector<std::vector<bool> > d_using_uniform_target_damping;
    std::vector<std::vector<double> > d_uniform_target_damping;

    /*
     * Anchor point information.
     */
    std::vector<std::vector<bool> > d_enable_anchor_points;

    struct AnchorSpec
    {
        bool is_anchor_point;
    };
    std::vector<std::vector<std::vector<AnchorSpec> > > d_anchor_spec_data;

    /*
     * Mass information for the pIB method.
     */
    std::vector<std::vector<bool> > d_enable_bdry_mass;

    struct BdryMassSpec
    {
        double bdry_mass, stiffness;
    };
    std::vector<std::vector<std::vector<BdryMassSpec> > > d_bdry_mass_spec_data;

    std::vector<std::vector<bool> > d_using_uniform_bdry_mass;
    std::vector<std::vector<double> > d_uniform_bdry_mass;

    std::vector<std::vector<bool> > d_using_uniform_bdry_mass_stiffness;
    std::vector<std::vector<double> > d_uniform_bdry_mass_stiffness;

    /*
     * Orthonormal directors for the generalized IB method.
     */
    std::vector<std::vector<std::vector<std::vector<double> > > > d_directors;

    /*
     * Instrumentation information.
     */
    std::vector<std::vector<bool> > d_enable_instrumentation;
    std::vector<std::vector<std::map<int, std::pair<int, int> > > > d_instrument_idx;

    /*
     * Source information.
     */
    std::vector<std::vector<bool> > d_enable_sources;
    std::vector<std::vector<std::map<int, int> > > d_source_idx;

    /*
     * Data required to specify connectivity information for visualization
     * purposes.
     */
    std::vector<unsigned int> d_global_index_offset;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBStandardInitializer
