// Filename: cRigidIBFEMethod.h
// Created on 14 Oct 2014 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of its
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

#ifndef included_cRigidIBFEMethod
#define included_cRigidIBFEMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stdbool.h>
#include <stddef.h>
#include <set>
#include <string>
#include <vector>

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "ibamr/cRigidIBStrategy.h"
#include "ibtk/libmesh_utilities.h"
#include "ibtk/FEDataManager.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace appu
{
template<int dim>
class VisItDataWriter;
}// namespace appu
}// namespace SAMRAI
namespace libMesh
{
class EquationSystems;
class Mesh;
class Point;
class System;
template <typename T>
class NumericVector;
template <typename T>
class PetscVector;
}// namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class cRigidIBFEMethod is a concrete cIBStrategy and RigidBodyStrategy
 * class which implements the motion of rigid bodies using the FE framework.
 */

class cRigidIBFEMethod : public cRigidIBStrategy
{

//////////////////////////////////////////////////////////////////////////////
	
public:
    static const std::string COORDS_SYSTEM_NAME;  
    static const std::string COORD_MAPPING_SYSTEM_NAME;
    static const std::string CONSTRAINT_FORCE_SYSTEM_NAME;
    static const std::string VELOCITY_SYSTEM_NAME;
    static const std::string CONSTRAINT_VELOCITY_SYSTEM_NAME;  
    static const std::string REGULATOR_SYSTEM_NAME;
      
    /*!
     * \brief Constructor.
     */
    cRigidIBFEMethod(
		const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        libMesh::Mesh* mesh,
        int max_level_number,
        bool register_for_restart = true);

    /*!
     * \brief Constructor.
     */
    cRigidIBFEMethod(
		const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::vector<libMesh::Mesh*>& meshes,
        int max_level_number,
        bool register_for_restart = true);
    
    /*!
     * \brief Destructor of the class.
     */
    ~cRigidIBFEMethod();
    
    /*!
     * Return a pointer to the finite element data manager object for the
     * specified part.
     */
    IBTK::FEDataManager* getFEDataManager(
        unsigned int part = 0) const;    
    
    /*!
     * Return the number of ghost cells required by the Lagrangian-Eulerian
     * interaction routines.
     */
    const SAMRAI::hier::IntVector<NDIM>& 
    getMinimumGhostCellWidth() const; 
    
    /*!
     * Typedef specifying interface for specifying constrained body velocities.
     */
    typedef void (*ConstrainedVelocityFcnPtr)(
        libMesh::NumericVector<double>& U_k,
	    const RigidDOFVector& U,
        libMesh::NumericVector<double>& X,
	    const Eigen::Vector3d& X_com,
        libMesh::EquationSystems* equation_systems,
        double data_time,
        void* ctx);
    
    typedef void (*ConstrainedCOMVelocityFcnPtr)(
        double data_time,
        Eigen::Vector3d& U_com,
	    Eigen::Vector3d& W_com);

    /*!
     * Struct encapsulating constrained velocity function data.
     */
    struct ConstrainedVelocityFcnData
    {
        ConstrainedVelocityFcnData(
	    ConstrainedVelocityFcnPtr    nodalvelfcn = NULL, 
	    ConstrainedCOMVelocityFcnPtr comvelfcn   = NULL,
	    void* ctx = NULL)
            : nodalvelfcn(nodalvelfcn), 
              comvelfcn(comvelfcn),
              ctx(ctx)
        {
        }

        ConstrainedVelocityFcnPtr    nodalvelfcn;
        ConstrainedCOMVelocityFcnPtr comvelfcn;	
        void* ctx;
    };

    /*!
     * Register a constrained body velocity function.
     */
    void 
    registerConstrainedVelocityFunction(
        ConstrainedVelocityFcnPtr    nodalvelfcn,
        ConstrainedCOMVelocityFcnPtr comvelfcn,
        void* ctx = NULL,
        unsigned int part = 0);

    /*!
     * Register a constrained body velocity function.
     */
    void 
    registerConstrainedVelocityFunction(
        const ConstrainedVelocityFcnData& data,
        unsigned int part = 0);
    
    /*!
     * Typedef specifying interface for coordinate mapping function.
     */
    typedef void (*CoordinateMappingFcnPtr)(
        libMesh::Point& X,
        const libMesh::Point& s,
        void* ctx);

    /*!
     * Struct encapsulating coordinate mapping function data.
     */
    struct CoordinateMappingFcnData
    {
        CoordinateMappingFcnData(CoordinateMappingFcnPtr fcn = NULL, void* ctx = NULL)
            : fcn(fcn), ctx(ctx)
        {
        }

        CoordinateMappingFcnPtr fcn;
        void* ctx;
    };

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
        CoordinateMappingFcnPtr fcn,
        void* ctx = NULL,
        unsigned int part = 0); 
	
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
        const CoordinateMappingFcnData& data,
        unsigned int part = 0);    
    
    /*!
     * Typedef specifying interface for Lagrangian body force distribution
     * function.
     */
    typedef IBTK::VectorMeshFcnPtr LagBodyForceFcnPtr;

    /*!
     * Struct encapsulating Lagrangian body force distribution data.
     */
    struct LagBodyForceFcnData
    {
        LagBodyForceFcnData(
            LagBodyForceFcnPtr fcn = NULL,
            const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
            void* ctx = NULL)
            : fcn(fcn), systems(systems), ctx(ctx)
        {
        }
        LagBodyForceFcnPtr fcn;
        std::vector<unsigned int> systems;
        void* ctx;
    };
    
    /*!
     * Register the (optional) function to compute body force distributions on
     * the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple body force functions
     * with this class.
     */
    void 
    registerLagBodyForceFunction(
        LagBodyForceFcnPtr fcn,
        const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
        void* ctx = NULL,
        unsigned int part = 0);

    /*!
     * Register the (optional) function to compute body force distributions on
     * the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple body force functions
     * with this class.
     */
    void 
    registerLagBodyForceFunction(
        const LagBodyForceFcnData& data, 
	    unsigned int part = 0);
    
    /*!
     * Typedef specifying interface for Lagrangian surface force distribution
     * function.
     */
    typedef IBTK::VectorSurfaceFcnPtr LagSurfaceForceFcnPtr;

    /*!
     * Struct encapsulating Lagrangian surface force distribution data.
     */
    struct LagSurfaceForceFcnData
    {
        LagSurfaceForceFcnData(
            LagSurfaceForceFcnPtr fcn = NULL,
            const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
            void* ctx = NULL)
            : fcn(fcn), systems(systems), ctx(ctx)
        {
        }

        LagSurfaceForceFcnPtr fcn;
        std::vector<unsigned int> systems;
        void* ctx;
    };

    /*!
     * Register the (optional) function to compute surface force distributions
     * on the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple surface force functions
     * with this class.
     */
    void 
    registerLagSurfaceForceFunction(
        LagSurfaceForceFcnPtr fcn,
        const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
        void* ctx = NULL,
        unsigned int part = 0);

    /*!
     * Register the (optional) function to compute surface force distributions
     * on the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple surface force functions
     * with this class.
     */
    void 
    registerLagSurfaceForceFunction(
        const LagSurfaceForceFcnData& data,
        unsigned int part = 0);    
    
    /*!
     * Setup the tag buffer.
     */
    virtual void 
    setupTagBuffer(
        SAMRAI::tbox::Array<int>& tag_buffer,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) const;

    /*!
     * \brief Method to prepare to advance data from current_time to new_time.
     */
    virtual void
    preprocessIntegrateData(
        double current_time,
        double new_time,
        int num_cycles); 
    
    /*!
     * \brief Method to clean up data following call(s) to integrateHierarchy().
     */
    virtual void
    postprocessIntegrateData(
        double current_time,
        double new_time,
        int num_cycles);
    
    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    virtual void 
    interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
            u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
            u_ghost_fill_scheds,
        double data_time);  
	
    /*!
     * \brief Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    virtual void
    eulerStep(
        double current_time,
        double new_time);

    /*!
     * \brief Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    virtual void
    midpointStep(
        double current_time,
        double new_time);  
    
    /*!
     * \brief Advance the positions of the Lagrangian structure using the 
     * trapezoidal rule.
     */
    virtual void
    trapezoidalStep(
        double current_time,
        double new_time);  
    
    /*!
     * \brief Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    virtual void 
    computeLagrangianForce(
        double data_time);
    
    /*!
     * \brief Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    virtual void
    spreadForce(
        int f_data_idx,
        IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
        f_prolongation_scheds,
        double data_time);
	
	// \{
	// The following methods are concrete implementation of cIBStrategy methods:
	//
	
	// \see cIBStrategy::spreadForce() methods.
	/*!
	 * \brief Spread the constraint Lagrangian force for all parts in the PetscMultiVec
	 * \f$L\f$ to the Cartesian grid at the specified time within the current
	 * time interval.
	 *
	 * \param scale Scales all the components of Lagrangian vector before spreading.
	 */
	virtual void
	spreadForce(
		int f_data_idx,
		Vec L,
		IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
		f_prolongation_scheds,
		double data_time,
		double scale = 1.0);
	
	/*!
	 * \brief Spread the constraint Lagrangian force for a specific part in the Vec
	 * \f$L\f$ to the Cartesian grid at the specified time within the current
	 * time interval.
	 *
	 * \param scale Scales the Lagrangian vector before spreading.
	 */
	virtual void
	spreadForce(
		int f_data_idx,
		const unsigned int part,
		Vec L,
		IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
		f_prolongation_scheds,
		double data_time,
		double scale = 1.0);
	
	/*!
	 * \brief Interpolate the Eulerian velocity to the curvilinear mesh for all parts
	 * in the PetscMultiVec \f$V\f$ at the specified time within the current
	 * time interval.
	 *
	 * \param scale Scales the Lagrangian vector after interpolating from Eulerian grid.
	 */
	virtual void
	interpolateVelocity(
	    int u_data_idx,
		Vec V,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
			u_synch_scheds,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
			u_ghost_fill_scheds,
		double data_time,
		double scale = 1.0);
	
	/*!
	 * \brief Interpolate the Eulerian velocity to the curvilinear mesh for a specific part
	 * in the Vec \f$V\f$ at the specified time within the current time interval.
	 *
	 ** \param scale Scales the Lagrangian vector after interpolating from Eulerian grid.
	 */
	virtual void
	interpolateVelocity(
		int u_data_idx,
		const unsigned int part,
		Vec V,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
			u_synch_scheds,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
			u_ghost_fill_scheds,
		double data_time,
		double scale = 1.0);
	// \}
	
	// \{
	// The following are the concrete implementation of RigidBodyStrategy methods:
	//
	
	// \see RigidBodyStrategy::computeNetRigidGeneralizedForce() methods.
	/*!
	 * \brief Compute total force and torque on the structure(s).
	 */
	virtual void
	computeNetRigidGeneralizedForce(
		const unsigned int part,
		Vec L,
		RigidDOFVector& F);
	
	// \see RigidBodyStrategy::setRigidBodyVelocity()
	/*!
	 * \brief Set the rigid body velocity at the nodal points
	 * contained in the Vec V.
	 *
	 */
	virtual void
	setRigidBodyVelocity(
		const unsigned int part,
		const RigidDOFVector& U,
		Vec V);
	
	// \see RigidBodyStrategy::getRigidBodyForce()
	/*!
	 * \brief Get the constraint rigid body force at the specified time within
	 * the current time interval.
	 */
	virtual void
	getRigidBodyForce(
		Vec* L,
		const double time);
	
	// \}
	
    /*!
     * Initialize FE data.  This method must be called prior to calling
     * IBHierarchyIntegrator::initializePatchHierarchy().
     */
    void initializeFEData();
    
    /*!
     * \brief Register Eulerian variables with the parent IBHierarchyIntegrator.
     */    
    virtual void
    registerEulerianVariables();
    
    /*!
     * \brief Register Eulerian refinement or coarsening algorithms with the parent
     * IBHierarchyIntegrator.
     */
    virtual void
    registerEulerianCommunicationAlgorithms();
    
    /*!
     * Initialize Lagrangian data corresponding to the given AMR patch hierarchy
     * at the start of a computation.  If the computation is begun from a
     * restart file, data may be read from the restart databases.
     *
     * A patch data descriptor is provided for the Eulerian velocity in case
     * initialization requires interpolating Eulerian data.  Ghost cells for
     * Eulerian data will be filled upon entry to this function.
     */
    void initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg,
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
            u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
            u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time);
    
    /*!
     * Register a load balancer and work load patch data index with the IB
     * strategy object.
     */
    virtual void 
    registerLoadBalancer(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
        int workload_data_idx); 

    /*!
     * Update work load estimates on each level of the patch hierarchy.
     */
    virtual void 
    updateWorkloadEstimates(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        int workload_data_idx);    
    
    /*!
     * Begin redistributing Lagrangian data prior to regridding the patch
     * hierarchy.
     */
    virtual void 
    beginDataRedistribution(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Complete redistributing Lagrangian data following regridding the patch
     * hierarchy.
     */
    virtual void 
    endDataRedistribution(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg); 
    
    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    virtual void initializeLevelData(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double init_data_time,
        bool can_be_refined,
        bool initial_time,
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
        bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    virtual void 
    resetHierarchyConfiguration(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int coarsest_level,
        int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to user-supplied feature detection criteria.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector
     */
    virtual void 
    applyGradientDetector(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double error_data_time,
        int tag_index,
        bool initial_time,
        bool uses_richardson_extrapolation_too);    
	
	/*!
	 * \brief Callbacks before INS is integrated.
	 */
	typedef void
	(*preprocessSolveFluidEqn_callbackfcn)(
	const double,
	const double,
	const int,
	void*);
	
    /*!
     * \brief Register any preprocess fluid solve callback functions.
     */
    void
    registerPreProcessSolveFluidEquationsCallBackFcn(
		preprocessSolveFluidEqn_callbackfcn callback,
        void* ctx);
    
    /*!
     * \brief Calculate any body forces for INS solver over here.
     */
    virtual void
    preprocessSolveFluidEquations(
        double current_time,
		double new_time,
		int cycle_num);
            
    /*!
     * Register a visIt data writer to output data files that
     * may be postprocessed with the visIt visualization tool.
     */
    void
    registerVisItDataWriter(
        SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);  
      
    /*!
     * \brief Get the level number on which structures reside.
     */
    int
    getStructuresLevelNumber();
	
    /*!
     * \brief Get the patch hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
    getPatchHierarchy();

   /*!
    *\brief Get the HierarchyMathOps object.
    */
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps>
    getHierarchyMathOps();
	
    /*!
     * \brief Override the putToDatabase method of the base Serializable class.
     */
    virtual void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
	
//////////////////////////////////////////////////////////////////////////////
protected:
   

//////////////////////////////////////////////////////////////////////////////    
private:
  
    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::vector<libMesh::Mesh*>& meshes,
        int max_level_number,
        bool register_for_restart);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, 
	    bool is_from_restart);

    /*!
     * \brief Get values from restart file.
     */
    void
    getFromRestart();
    
    /*!
     * \brief Compute dX = X - s, useful mainly for visualization purposes.
     */
    void 
    updateCoordinateMapping(
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
     * \brief Compute center of mass and moment of inertia of structures.
     */
    void
    computeCOMandMOIOfStructures(
        std::vector<Eigen::Vector3d>& center_of_mass,
        std::vector<Eigen::Matrix3d>& moment_of_inertia,
        std::vector<libMesh::PetscVector<double>*> X);

    /*!
     * \brief Fill the rotation matrix.
     * \param rot_mat Matrix to set. 
     * \param dt Time interval of rotation. 
     */
    void
    setRotationMatrix(
        const std::vector<Eigen::Vector3d>& rot_vel,
        std::vector<Eigen::Matrix3d>& rot_mat,
        const double dt);
    
    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log;

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;
    bool d_is_initialized;
	
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
    
    /*
     * The current time step interval.
     */
    double d_current_time, d_new_time, d_half_time;
    
    /*
     * Number of nodes of rigid structures.
     */
    std::vector<unsigned int> d_num_nodes;

    /*
     * FE data associated with this object.
     */
    std::vector<libMesh::Mesh*> d_meshes;
    std::vector<libMesh::EquationSystems*> d_equation_systems;

    std::vector<IBTK::FEDataManager*> d_fe_data_managers;
    SAMRAI::hier::IntVector<NDIM> d_ghosts;
    std::vector<libMesh::System*> d_X_systems, d_U_systems, d_L_systems, d_U_constrained_systems;
    std::vector<libMesh::PetscVector<double>*> d_X_current_vecs, d_X_new_vecs, d_X_half_vecs,
        d_X_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_current_vecs, d_U_new_vecs, d_U_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_L_current_vecs, d_L_new_vecs, d_L_half_vecs,
		d_L_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_constrained_current_vecs, d_U_constrained_new_vecs,
        d_U_constrained_half_vecs;
    bool d_fe_data_initialized;
	
	/*!
	 * PETSc wrappers for rigid body force.
	 */
	std::vector<Vec> d_vL_new, d_vL_current;
	Vec d_mv_L_new, d_mv_L_current;
	
    /*
     * Method paramters.
     */
    IBTK::FEDataManager::InterpSpec d_interp_spec;
    IBTK::FEDataManager::SpreadSpec d_spread_spec;
    libMeshEnums::FEFamily d_fe_family;
    libMeshEnums::Order d_fe_order;
    libMeshEnums::QuadratureType d_quad_type;
    libMeshEnums::Order d_quad_order;
    bool d_use_consistent_mass_matrix;
     
    /*
     * Nonuniform load balancing data structures.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;
    int d_workload_idx;
    
    /*
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<CoordinateMappingFcnData> d_coordinate_mapping_fcn_data;
    
    /*
     * Functions used to compute additional body and surface forces on the
     * Lagrangian mesh.
     */
    std::vector<LagBodyForceFcnData> d_lag_body_force_fcn_data;
    std::vector<LagSurfaceForceFcnData> d_lag_surface_force_fcn_data;   
    
    /*
     * Collection of all systems required to evaluate various quantities.
     */
    std::vector<std::set<unsigned int> > d_fcn_systems, d_body_fcn_systems,
        d_surface_fcn_systems;    

    /*!
     * Functions to set constrained velocities of the structures.
     */
    std::vector<ConstrainedVelocityFcnData> d_constrained_velocity_fcn_data;
	
    /*!
     * Pre and post fluid solve call back functions and contexts.
     */
    std::vector<preprocessSolveFluidEqn_callbackfcn> d_prefluidsolve_callback_fcns;
    std::vector<void*> d_prefluidsolve_callback_fcns_ctx;
    
    /*!
     * Eulerian variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_eul_lambda_var;
    int d_eul_lambda_idx;
    
    /*!
     * The object used to write out data for postprocessing by the visIt
     * visualization tool.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;

    /*!
     * Control printing of S[lambda]
     */
    bool d_output_eul_lambda;
 
};//cRigidIBFEMethod
}// namespace IBAMR


//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_cRigidIBFEMethod
