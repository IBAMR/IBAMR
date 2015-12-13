// Filename: CIBMethod.h
// Created on 21 Apr 2015 by Amneet Bhalla and Bakytzhan Kallemov 
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith.
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

#ifndef included_CIBMethod
#define included_CIBMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "RobinBcCoefStrategy.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/IBMethod.h"
#include "ibamr/IBMethod.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"

namespace IBTK
{
class HierarchyMathsOps;
} // namespace IBTK
namespace IBAMR
{
class CIBStandardInitializer;
class CIBStaggeredStokesOperator;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CIBFEMethod is a concrete CIBStrategy and IBMethod
 * class which implements the motion of rigid bodies using the constraint
 * formulation. The immersed structure is discretized using standard IB
 * markers.
 */

class CIBMethod : public IBAMR::IBMethod, public IBAMR::CIBStrategy
{

public:
    /*!
     * \brief Constructor of the class.
     */
    CIBMethod(const std::string& object_name,
              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
              const int no_structures = 1,
              bool register_for_restart = true);

    /*!
     * \brief Destructor of the class.
     */
    ~CIBMethod();

    /*!
     * \brief Typedef specifying interface for specifying constrained body velocities.
     */
    typedef void (*ConstrainedNodalVelocityFcnPtr)(const unsigned part,
						   Vec U_k,
                                                   const RigidDOFVector& U,
                                                   Vec X,
                                                   const Eigen::Vector3d& X_com,
                                                   double data_time,
                                                   void* ctx,
						   IBAMR::CIBMethod* cib_method);

    typedef void (*ConstrainedCOMVelocityFcnPtr)(double data_time, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com);

    /*!
     * \brief Typedef specifying interface for specifying net external force and torque on structures.
     */
    typedef void (*ExternalForceTorqueFcnPtr)(double data_time, Eigen::Vector3d& F, Eigen::Vector3d& T);

    /*!
     * \brief Callbacks before INS is integrated.
     */
    typedef void (*preprocessSolveFluidEqn_callbackfcn)(const double, const double, const int, void*);

    /*!
     * \brief Typedef specifying interface for specifying slip velocities.
     */
    typedef void (*VelocityDeformationFunctionPtr)(Vec V,
                                                   Vec X,
                                                   const std::vector<Eigen::Vector3d>& X_com,
						   IBAMR::CIBMethod* d_cib_method);

    /*!
     * \brief Struct encapsulating constrained velocity functions data.
     */
    struct ConstrainedVelocityFcnsData
    {
        ConstrainedVelocityFcnsData(ConstrainedNodalVelocityFcnPtr nodalvelfcn = NULL,
                                    ConstrainedCOMVelocityFcnPtr comvelfcn = NULL,
                                    void* ctx = NULL)
            : nodalvelfcn(nodalvelfcn), comvelfcn(comvelfcn), ctx(ctx)
        {
            // intentionally left blank
        }
        ConstrainedNodalVelocityFcnPtr nodalvelfcn;
        ConstrainedCOMVelocityFcnPtr comvelfcn;
        void* ctx;
    };

    /*!
     * \brief Struct encapsulating external force and torque function data.
     */
    struct ExternalForceTorqueFcnData
    {
        ExternalForceTorqueFcnData(ExternalForceTorqueFcnPtr forcetorquefcn = NULL, void* ctx = NULL)
            : forcetorquefcn(forcetorquefcn), ctx(ctx)
        {
            // intentionally left blank
        }
        ExternalForceTorqueFcnPtr forcetorquefcn;
        void* ctx;
    };

    /*!
     * \brief Register a constrained body velocity function.
     */
    void registerConstrainedVelocityFunction(ConstrainedNodalVelocityFcnPtr nodalvelfcn,
                                             ConstrainedCOMVelocityFcnPtr comvelfcn,
                                             void* ctx = NULL,
                                             unsigned int part = 0);
    /*!
     * \brief Register a constrained body velocity function data.
     */
    void registerConstrainedVelocityFunction(const ConstrainedVelocityFcnsData& data, unsigned int part = 0);

    /*!
     * \brief register deformation velocity.
     */
    void registerRigidBodyVelocityDeformationFunction(VelocityDeformationFunctionPtr VelDefFun);

    /*!
     * \brief Register an external force and torque function.
     */
    void registerExternalForceTorqueFunction(ExternalForceTorqueFcnPtr forcetorquefcn,
                                             void* ctx = NULL,
                                             unsigned int part = 0);
    /*!
     * \brief Register an external force and torque function data.
     */
    void registerExternalForceTorqueFunction(const ExternalForceTorqueFcnData& data, unsigned int part = 0);

    /*!
     * \brief Get the level on which the structures reside.
     */
    int getStructuresLevelNumber() const;

    /*!
     * \brief Get the structure handle to which this Lagrangian index belongs.
     */
    int getStructureHandle(const int lag_idx) const;

    /*!
     * \brief Register any preprocess fluid solve callback functions.
     */
    void registerPreProcessSolveFluidEquationsCallBackFcn(preprocessSolveFluidEqn_callbackfcn callback, void* ctx);

    /*!
     * \brief Calculate any body forces for INS solver over here.
     */
    void preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num);

    /*!
     * \brief Register Eulerian variables with the parent IBHierarchyIntegrator.
     */
    void registerEulerianVariables();

    /*!
     * \brief Register Eulerian refinement or coarsening algorithms with the parent
     * IBHierarchyIntegrator.
     */
    void registerEulerianCommunicationAlgorithms();

    /*!
     * \brief Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * \brief Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * \brief Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                             int level_number,
                             double init_data_time,
                             bool can_be_refined,
                             bool initial_time,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                             bool allocate_data);

    /*!
     * \brief Initialize Lagrangian data corresponding to the given AMR patch hierarchy
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
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time);

    /*!
    * \brief Interpolate the Eulerian velocity to the curvilinear mesh at the
    * specified time within the current time interval.
    */
    void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time);

    /*!
     * \brief Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time);

    /*!
     * \brief Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void eulerStep(double current_time, double new_time);

    /*!
     * \brief Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void midpointStep(double current_time, double new_time);

    /*!
     * \brief Advance the positions of the Lagrangian structure using the
     * trapezoidal rule.
     */
    void trapezoidalStep(double current_time, double new_time);

    /*!
     * \brief Register VisIt data writer to output data files that
     * may be postprocessed with the VisIt visualization tool.
     */
    void registerVisItDataWriter(SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * \brief Override the putToDatabase method of the base Serializable class.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    // \{
    // The following are the concrete implementation of CIBStrategy methods:
    //

    // \see CIBStrategy::setConstraintForce() method.
    /*!
     * \brief Set the constraint force in the internal data structures of the
     * class.
     */
    void setConstraintForce(Vec L, const double data_time, const double scale = 1.0);

    // \see CIBStrategy::getConstraintForce()
    /*!
     * \brief Get the constraint rigid body force at the specified time within
     * the current time interval.
     */
    void getConstraintForce(Vec* L, const double data_time);

    // \see CIBStrategy::getFreeRigidVelocities()
    /*!
     * \brief Get the free rigid velocities (DOFs) at the specified time within
     * the current time interval.
     */
    virtual void getFreeRigidVelocities(Vec* U, const double data_time);

    // \see CIBStrategy::getNetExternalForceTorque()
    /*!
     * \brief Get net external force and torque at the specified time within
     * the current time interval.
     */
    virtual void getNetExternalForceTorque(Vec* F, const double data_time);

    // \see CIBStrategy::subtractMeanConstraintForce()
    /*!
     * \brief Subtract the mean of constraint force from the background Eulerian
     * grid.
     */
    void subtractMeanConstraintForce(Vec L, int f_data_idx, const double scale = 1.0);

    // \see CIBStrategy::setInterpolatedVelocityVector() method
    /*!
     * \brief Prepare the CIBMethod class to get the interpolated fluid
     * velocity.
     */
    void setInterpolatedVelocityVector(Vec V, const double data_time);

    // \see CIBStrategy::setInterpolatedVelocityVector() method
    /*!
     * \brief Get interpolated velocity from the Eulerian grid.
     */
    void getInterpolatedVelocity(Vec V, const double data_time, const double scale = 1.0);

    // \see CIBStrategy::computeMobilityRegularization method
    /*!
     * \brief Compute regularization vector for the mobility problem.
     *
     */
    void computeMobilityRegularization(Vec D, Vec L, const double scale = 1.0);

    // \see CIBStrategy::getNumberOfNodes method
    /*!
     * \brief Get number of nodes for a particular structure registered with CIBMethod.
     */
    unsigned int getNumberOfNodes(const unsigned int part) const;

    // \see CIBStrategy::setRigidBodyVelocity method
    /*!
     * \brief Set the rigid body velocity at the nodal points
     * contained in the Vec V.
     */
    void setRigidBodyVelocity(Vec U, Vec V, const std::vector<bool>& skip_comp, const bool isHalfTimeStep=true);

    // \see CIBStrategy::setRigidBodyDeformationVelocity method.
    /*!
     * \brief set deformation velocity.
     */
    void setRigidBodyDeformationVelocity(Vec W);

    //void setRigidBodyDeformationVelocity(Vec W);


    // \see CIBStrategy::computeNetRigidGeneralizedForce() method.
    /*!
     * \brief Compute total force and torque on the rigid structure(s).
     */
    void computeNetRigidGeneralizedForce(Vec L, Vec F, const std::vector<bool>& skip_comp, const bool isHalfTimeStep=true);


    void copyAllArrayToVec(Vec b,
			   double* array,
			   const std::vector<unsigned>& all_rhs_struct_ids,
			   const int data_depth);
    
    void copyAllVecToArray(Vec b,
			   double* array,
			   const std::vector<unsigned>& all_rhs_struct_ids,
			   const int data_depth);

    // \see CIBStrategy::constructMobilityMatrix() method.
    /*!
     * \brief Generate dense mobility matrix for the prototypical structures
     * identified by their indices.
     */
    void constructMobilityMatrix(std::map<std::string, double*>& managed_mat_map,
				 std::map<std::string, MobilityMatrixType>& managed_mat_type_map,
				 std::map<std::string, std::vector<unsigned> >& managed_mat_prototype_id_map,
				 std::map<std::string, unsigned int>& managed_mat_nodes_map,
				 std::map<std::string, std::pair<double, double> >& managed_mat_scale_map,
				 std::map<std::string, int>& managed_mat_proc_map,
				 const double* grid_dx,
				 const double* domain_extents,
				 const bool initial_time,
				 double rho,
				 double mu,
				 double f_periodic_corr);

    // \see CIBStrategy::constructBodyMobilityMatrix() method.
    /*!
     * \brief Generate dense body mobility matrix for the free structures
     * identified by their ids.
     */
    void constructKinematicMatrix(double* kinematic_mat,
				  const std::vector<unsigned>& prototype_struct_ids,
				  const bool initial_time,
				  const int managing_rank);
    
   //\see CIBStrategy:: rotateArrayInitalBodyFrame method.
    /*!
     * \brief rotate vector using stored structures rotation matrix to/from the reference frame of the structures at initial time.
     *
     */
    void rotateArrayInitalBodyFrame(double* array, 
				    const std::vector<unsigned>& struct_ids,
				    const bool isTranspose,
				    const int managing_rank,
				    const bool BodyMobility = false);


    /*!
     * \brief register CIBStandartInitializer 
     */
    void registerStandardInitializer(SAMRAI::tbox::Pointer<IBAMR::CIBStandardInitializer> ib_initializer);
    
    /*!
     * \brief returns IBStandartInitializer 
     */
    SAMRAI::tbox::Pointer<IBAMR::CIBStandardInitializer> getStandardInitializer();
    
    /*!
     * \brief Compute initial center of massof structures.
     */
    void computeInitialCOMOfStructures(std::vector<Eigen::Vector3d>& center_of_mass);

    /*!
     * \brief Compute initial moment of inertia of structures.
     */
    void computeInitialMOIOfStructures(std::vector<Eigen::Matrix3d>& moment_of_inertia);

    /*!
     * Function to determine whether regridding should occur at the current time
     * step.
     */
    bool flagRegrid() const;

    /*
     * \brief Function to fill ghost cells using values from the interior domain
     * and boundary conditions.
     */
    void fillGhostCells(int in, const double time);

    /*
     * Set velocity boundary conditions
     */
    void setVelocityBC(std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> *u_bc_coefs);

    /*
     * Set velocity physical boundary options
     */
    void setVelocityPhysBdryOp(IBTK::RobinPhysBdryPatchStrategy* u_phys_bdry_op);

    //////////////////////////////////////////////////////////////////////////////

protected:
    //////////////////////////////////////////////////////////////////////////////

private:
    /*!
     * \brief Set additional values from input database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Get values from restart file.
     */
    void getFromRestart();

    /*!
     * \brief Set regularization weight for Lagrangian markers.
     */
    void setRegularizationWeight(const int level_number);

    /*!
     * \brief Set initial Lambda for Lagrangian markers.
     */
    void setInitialLambda(const int level_number);

    /*!
     * \brief pointer to CIBStandartInitializer.
     */
    SAMRAI::tbox::Pointer<IBAMR::CIBStandardInitializer> d_ib_initializer;

    /*!
     * Functions to set constrained velocities of the structures.
     */
    std::vector<ConstrainedVelocityFcnsData> d_constrained_velocity_fcns_data;

    /*!
     * Functions to set net external force and torque on free moving structures.
     */
    std::vector<ExternalForceTorqueFcnData> d_ext_force_torque_fcn_data;

    /*!
     * Pre and post fluid solve call back functions and contexts.
     */
    std::vector<preprocessSolveFluidEqn_callbackfcn> d_prefluidsolve_callback_fcns;
    std::vector<void*> d_prefluidsolve_callback_fcns_ctx;

    /*!
     * Booleans to control spreading constraint force and interpolating
     * to Lagrangian velocities.
     */
    bool d_constraint_force_is_initialized, d_lag_velvec_is_initialized;

    /*!
     * Boolean to flag if time integrator needs regriding
     */
    bool d_time_integrator_needs_regrid;

    /*!
     * Eulerian variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_eul_lambda_var;
    int d_eul_lambda_idx;

    /*!
     * Vector of Lagrnagian indices of all structures.
     */
    std::vector<std::pair<int, int> > d_struct_lag_idx_range;

    /*!
     * PETSc wrappers for free rigid body velocities.
     */
    Vec d_mv_U, d_mv_F;

    /*!
     * The object used to write out data for postprocessing by the visIt
     * visualization tool.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;

    /*!
     * Control printing of S[lambda]
     */
    bool d_output_eul_lambda, d_output_all_lambdas;

    /*!
     * Printing Lagrange multiplier.
     */
    int d_lambda_dump_interval;
    std::ofstream d_lambda_stream;
    std::ofstream d_netlambda_stream;

    std::vector<std::string> d_reg_filename;
    std::vector<std::string> d_lambda_filename;
    VelocityDeformationFunctionPtr d_VelDefFun;

    // velocity boundary coefficients
    // vector<RobinBcCoefStrategy<NDIM>*> *d_u_bc_coefs;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> *d_u_bc_coefs;
    IBTK::RobinPhysBdryPatchStrategy* d_u_phys_bdry_op;   

}; // CIBMethod
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CIBMethod
