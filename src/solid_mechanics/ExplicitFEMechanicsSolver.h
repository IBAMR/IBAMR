// Filename: ExplicitFEMechanicsSolver.h
// Created on 12 Mar 2012 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_ExplicitFEMechanicsSolver
#define included_ExplicitFEMechanicsSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

// LIBMESH INCLUDES
#define LIBMESH_REQUIRE_SEPARATE_NAMESPACE
#include <../base/variable.h>
#include <enum_order.h>
#include <enum_quadrature_type.h>
#include <equation_systems.h>
#include <linear_solver.h>
#include <mesh.h>
#include <petsc_vector.h>
#include <sparse_matrix.h>

// PETSC INCLUDES
#include <petscsys.h>

// SAMRAI INCLUDES
#include <tbox/Serializable.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class ExplicitFEMechanicsSolver is a simple explicit mechanics solver
 * using a purely displacement-based formulation.  Multi-part models are
 * supported, but we do not attempt to handle contact between structures.
 */
class ExplicitFEMechanicsSolver
    : public SAMRAI::tbox::Serializable

{
public:
    static const std::string        COORDS_SYSTEM_NAME;
    static const std::string COORD_MAPPING_SYSTEM_NAME;
    static const std::string         FORCE_SYSTEM_NAME;
    static const std::string      VELOCITY_SYSTEM_NAME;
    static const std::string     F_DIL_BAR_SYSTEM_NAME;

    /*!
     * \brief The libMesh boundary IDs to use for specifying essential boundary
     * conditions.
     *
     * \todo Move these to a common header file.
     */
//  static const short int     NORMAL_DIRICHLET_BDRY_ID = 256;
//  static const short int TANGENTIAL_DIRICHLET_BDRY_ID = 512;
    static const short int            DIRICHLET_BDRY_ID = 256 | 512;

    /*!
     * \brief Constructor.
     */
    ExplicitFEMechanicsSolver(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        libMesh::Mesh* mesh,
        bool register_for_restart=true);

    /*!
     * \brief Constructor.
     */
    ExplicitFEMechanicsSolver(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::vector<libMesh::Mesh*>& meshes,
        bool register_for_restart=true);

    /*!
     * \brief Destructor.
     */
    ~ExplicitFEMechanicsSolver();

    /*!
     * Return a pointer to the EquationSystems object for the specified part.
     */
    libMesh::EquationSystems*
    getEquationSystems(
        unsigned int part=0) const;

    /*!
     * Typedef specifying interface for coordinate mapping function.
     */
    typedef
    void
    (*CoordinateMappingFcnPtr)(
        libMesh::Point& X,
        const libMesh::Point& s,
        void* ctx);

    /*!
     * Register the (optional) function used to initialize the physical
     * coordinates from the Lagrangian coordinates.
     *
     * \note If no function is provided, the initial physical coordinates are
     * taken to be the same as the Lagrangian coordinate system, i.e., the
     * initial coordinate mapping is assumed to be the identity mapping.
     */
    void
    registerInitialCoordinateMappingFunction(
        CoordinateMappingFcnPtr coordinate_mapping_fcn,
        void* coordinate_mapping_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for PK1 stress tensor function.
     */
    typedef
    void
    (*PK1StressFcnPtr)(
        libMesh::TensorValue<double>& PP,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        libMesh::NumericVector<double>& X_vec,
        const std::vector<libMesh::NumericVector<double>*>& system_data,
        double time,
        void* ctx);

    /*!
     * Register the (optional) function to compute the first Piola-Kirchhoff
     * stress tensor, used to compute the forces on the Lagrangian finite
     * element mesh.
     */
    void
    registerPK1StressTensorFunction(
        PK1StressFcnPtr PK1_stress_fcn,
        std::vector<unsigned int> PK1_stress_fcn_systems=std::vector<unsigned int>(),
        void* PK1_stress_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for Lagrangian body force distribution
     * function.
     */
    typedef
    void
    (*LagBodyForceFcnPtr)(
        libMesh::VectorValue<double>& F,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        libMesh::NumericVector<double>& X_vec,
        const std::vector<libMesh::NumericVector<double>*>& system_data,
        double time,
        void* ctx);

    /*!
     * Register the (optional) function to compute body force distributions on
     * the Lagrangian finite element mesh.
     */
    void
    registerLagBodyForceFunction(
        LagBodyForceFcnPtr lag_body_force_fcn,
        std::vector<unsigned int> lag_body_force_fcn_systems=std::vector<unsigned int>(),
        void* lag_body_force_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for Lagrangian pressure force distribution
     * function.
     */
    typedef
    void
    (*LagPressureFcnPtr)(
        double& P,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        unsigned short int side,
        libMesh::NumericVector<double>& X_vec,
        const std::vector<libMesh::NumericVector<double>*>& system_data,
        double time,
        void* ctx);

    /*!
     * Register the (optional) function to compute surface pressure
     * distributions on the Lagrangian finite element mesh.
     */
    void
    registerLagPressureFunction(
        LagPressureFcnPtr lag_pressure_fcn,
        std::vector<unsigned int> lag_pressure_fcn_systems=std::vector<unsigned int>(),
        void* lag_pressure_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for Lagrangian surface force distribution
     * function.
     */
    typedef
    void
    (*LagSurfaceForceFcnPtr)(
        libMesh::VectorValue<double>& F,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        unsigned short int side,
        libMesh::NumericVector<double>& X_vec,
        const std::vector<libMesh::NumericVector<double>*>& system_data,
        double time,
        void* ctx);

    /*!
     * Register the (optional) function to compute surface force distributions
     * on the Lagrangian finite element mesh.
     */
    void
    registerLagSurfaceForceFunction(
        LagSurfaceForceFcnPtr lag_surface_force_fcn,
        std::vector<unsigned int> lag_surface_force_fcn_systems=std::vector<unsigned int>(),
        void* lag_surface_force_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void
    preprocessIntegrateData(
        double current_time,
        double new_time);

    /*!
     * Method to advance data from current_time to new_time.
     */
    void
    integrateData(
        double current_time,
        double new_time);

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateData(
        double current_time,
        double new_time);

    /*!
     * Initialize FE data.
     */
    void
    initializeFEData();

    /*!
     * Write out object state to the given database.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*
     * \brief Compute the projected dilatational strain F_dil_bar.
     */
    void
    computeProjectedDilatationalStrain(
        libMesh::NumericVector<double>& F_dil_bar_vec,
        libMesh::NumericVector<double>& X_vec,
        unsigned int part);

    /*
     * \brief Compute the interior elastic density, possibly splitting off the
     * normal component of the transmission force along the physical boundary of
     * the Lagrangian structure.
     */
    void
    computeInteriorForceDensity(
        libMesh::NumericVector<double>& G_vec,
        libMesh::NumericVector<double>& X_vec,
        libMesh::NumericVector<double>* F_dil_bar_vec,
        double time,
        unsigned int part);

    /*!
     * \brief Initialize the physical coordinates using the supplied coordinate
     * mapping function.  If no function is provided, the initial coordinates
     * are taken to be the Lagrangian coordinates.
     */
    void
    initializeCoordinates(
        unsigned int part);

    /*!
     * \brief Compute dX = X - s, useful mainly for visualization purposes.
     */
    void
    updateCoordinateMapping(
        unsigned int part);

    /*!
     * \return Pointers to a linear solver and sparse matrix corresponding to a
     * L2 projection operator.
     */
    std::pair<libMesh::LinearSolver<double>*,libMesh::SparseMatrix<double>*>
    buildL2ProjectionSolver(
        const std::string& system_name,
        unsigned int part=0,
        libMeshEnums::QuadratureType quad_type=QGAUSS,
        libMeshEnums::Order quad_order=FIFTH);

    /*!
     * \return Pointer to vector representation of diagonal L2 mass matrix.
     */
    libMesh::NumericVector<double>*
    buildDiagonalL2MassMatrix(
        const std::string& system_name,
        unsigned int part=0);

    /*!
     * \brief Set U to be the L2 projection of F.
     */
    bool
    computeL2Projection(
        libMesh::NumericVector<double>& U,
        libMesh::NumericVector<double>& F,
        const std::string& system_name,
        unsigned int part=0,
        bool consistent_mass_matrix=true,
        libMeshEnums::QuadratureType quad_type=QGAUSS,
        libMeshEnums::Order quad_order=FIFTH,
        double tol=1.0e-6,
        unsigned int max_its=100);

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log;

    /*
     * Indicates whether the FE data have been initialized.
     */
    bool d_is_initialized;

    /*
     * The current time step interval.
     */
    double d_current_time, d_new_time, d_half_time;

    /*
     * FE data associated with this object.
     */
    std::vector<libMesh::Mesh*> d_meshes;
    std::vector<libMesh::EquationSystems*> d_equation_systems;

    const unsigned int d_num_parts;
    std::vector<libMesh::System*> d_X_systems, d_U_systems, d_F_systems, d_F_dil_bar_systems;
    std::vector<libMesh::PetscVector<double>*> d_X_current_vecs, d_X_new_vecs, d_X_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_current_vecs, d_U_new_vecs, d_U_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_F_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_F_dil_bar_half_vecs;

    /*
     * Method paramters.
     */
    bool d_use_consistent_mass_matrix;
    bool d_use_Fbar_projection;
    libMeshEnums::FEFamily d_fe_family;
    libMeshEnums::Order d_fe_order;
    libMeshEnums::FEFamily d_F_dil_bar_fe_family;
    libMeshEnums::Order d_F_dil_bar_fe_order;
    libMeshEnums::QuadratureType d_quad_type;
    libMeshEnums::Order d_quad_order;

    /*
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<CoordinateMappingFcnPtr> d_coordinate_mapping_fcns;
    std::vector<void*> d_coordinate_mapping_fcn_ctxs;

    /*
     * Functions used to compute the first Piola-Kirchhoff stress tensor.
     */
    std::vector<PK1StressFcnPtr> d_PK1_stress_fcns;
    std::vector<std::vector<unsigned int> > d_PK1_stress_fcn_systems;
    std::vector<void*> d_PK1_stress_fcn_ctxs;

    /*
     * Functions used to compute additional body and surface forces on the
     * Lagrangian mesh.
     */
    std::vector<LagBodyForceFcnPtr> d_lag_body_force_fcns;
    std::vector<std::vector<unsigned int> > d_lag_body_force_fcn_systems;
    std::vector<void*> d_lag_body_force_fcn_ctxs;

    std::vector<LagPressureFcnPtr> d_lag_pressure_fcns;
    std::vector<std::vector<unsigned int> > d_lag_pressure_fcn_systems;
    std::vector<void*> d_lag_pressure_fcn_ctxs;

    std::vector<LagSurfaceForceFcnPtr> d_lag_surface_force_fcns;
    std::vector<std::vector<unsigned int> > d_lag_surface_force_fcn_systems;
    std::vector<void*> d_lag_surface_force_fcn_ctxs;

    /*
     * The (uniform) mass density of the structure in the reference
     * configuration.
     */
    double d_rho0;

    /*
     * Linear solvers and related data.
     */
    std::vector<std::map<std::string,libMesh::LinearSolver<double>*> > d_L2_proj_solver;
    std::vector<std::map<std::string,libMesh::SparseMatrix<double>*> > d_L2_proj_matrix;
    std::vector<std::map<std::string,libMesh::NumericVector<double>*> > d_L2_proj_matrix_diag;
    std::vector<std::map<std::string,libMeshEnums::QuadratureType> > d_L2_proj_quad_type;
    std::vector<std::map<std::string,libMeshEnums::Order> > d_L2_proj_quad_order;

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * A boolean value indicating whether the class is registered with the
     * restart database.
     */
    bool d_registered_for_restart;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    ExplicitFEMechanicsSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    ExplicitFEMechanicsSolver(
        const ExplicitFEMechanicsSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    ExplicitFEMechanicsSolver&
    operator=(
        const ExplicitFEMechanicsSolver& that);

    /*!
     * Implementation of class constructor.
     */
    void
    commonConstructor(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::vector<libMesh::Mesh*>& meshes,
        bool register_for_restart);

    /*!
     * Read input values from a given database.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void
    getFromRestart();
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/ExplicitFEMechanicsSolver.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ExplicitFEMechanicsSolver
