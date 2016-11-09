// Filename: CIBFEMethod.h
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

#ifndef included_IBAMR_CIBFEMethod
#define included_IBAMR_CIBFEMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <set>
#include <stdbool.h>
#include <stddef.h>
#include <string>
#include <vector>

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/IBFEMethod.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class HierarchyMathsOps;
}
namespace SAMRAI
{
namespace appu
{
template <int dim>
class VisItDataWriter;
} // namespace appu
} // namespace SAMRAI
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
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CIBFEMethod is a concrete CIBStrategy and IBFEMethod
 * class which implements the motion of rigid bodies using the constraint
 * formulation. The immersed structure is discretized using finite-element
 * mesh.
 */

class CIBFEMethod : public IBFEMethod, public CIBStrategy
{
    //////////////////////////////////////////////////////////////////////////////

public:
    static const std::string CONSTRAINT_VELOCITY_SYSTEM_NAME;

    /*!
     * \brief Constructor.
     */
    CIBFEMethod(const std::string& object_name,
                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                libMesh::Mesh* mesh,
                int max_level_number,
                bool register_for_restart = true);

    /*!
     * \brief Constructor.
     */
    CIBFEMethod(const std::string& object_name,
                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                const std::vector<libMesh::Mesh*>& meshes,
                int max_level_number,
                bool register_for_restart = true);

    /*!
     * \brief Destructor of the class.
     */
    ~CIBFEMethod();

    /*!
     * Typedef specifying interface for specifying constrained body velocities.
     */
    typedef void (*ConstrainedNodalVelocityFcnPtr)(libMesh::NumericVector<double>& U_k,
                                                   const RigidDOFVector& U,
                                                   libMesh::NumericVector<double>& X,
                                                   const Eigen::Vector3d& X_com,
                                                   libMesh::EquationSystems* equation_systems,
                                                   double data_time,
                                                   void* ctx);

    typedef void (*ConstrainedCOMVelocityFcnPtr)(double data_time, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com);

    /*!
     * \brief Typedef specifying interface for specifying net external force and torque on structures.
     */
    typedef void (*ExternalForceTorqueFcnPtr)(double data_time, Eigen::Vector3d& F, Eigen::Vector3d& T);

    /*!
     * Struct encapsulating constrained velocity functions data.
     */
    struct ConstrainedVelocityFcnsData
    {
        ConstrainedVelocityFcnsData(ConstrainedNodalVelocityFcnPtr nodalvelfcn = NULL,
                                    ConstrainedCOMVelocityFcnPtr comvelfcn = NULL,
                                    void* ctx = NULL)
            : nodalvelfcn(nodalvelfcn), comvelfcn(comvelfcn), ctx(ctx)
        {
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
     * Register a constrained body velocity function.
     */
    void registerConstrainedVelocityFunction(ConstrainedNodalVelocityFcnPtr nodalvelfcn,
                                             ConstrainedCOMVelocityFcnPtr comvelfcn,
                                             void* ctx = NULL,
                                             unsigned int part = 0);

    /*!
     * Register a constrained body velocity function data.
     */
    void registerConstrainedVelocityFunction(const ConstrainedVelocityFcnsData& data, unsigned int part = 0);

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
     * \brief Method to prepare to advance data from current_time to new_time.
     */
    virtual void preprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * \brief Method to clean up data following call(s) to integrateHierarchy().
     */
    virtual void postprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * \brief Indicate whether L2-projection is to be performed for velocity
     * interpolation.
     */
    bool setComputeVelL2Projection(const bool compute_L2_projection);

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    virtual void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time);

    /*!
     * \brief Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    virtual void eulerStep(double current_time, double new_time);

    /*!
     * \brief Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    virtual void midpointStep(double current_time, double new_time);

    /*!
     * \brief Advance the positions of the Lagrangian structure using the
     * trapezoidal rule.
     */
    virtual void trapezoidalStep(double current_time, double new_time);

    /*!
     * \brief Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    virtual void computeLagrangianForce(double data_time);

    /*!
     * \brief Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    virtual void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time);

    // \{
    // The following are the concrete implementation of CIBStrategy methods:
    //

    // \see CIBStrategy::setConstraintForce() method.
    /*!
     * \brief Set the constraint force in the internal data structures of the
     * class.
     */
    virtual void setConstraintForce(Vec L, const double data_time, const double scale = 1.0);

    // \see CIBStrategy::getConstraintForce()
    /*!
     * \brief Get the constraint rigid body force at the specified time within
     * the current time interval.
     */
    virtual void getConstraintForce(Vec* L, const double data_time);

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

    // \see CIBStrategy::subtractMeanConstraintForce().
    /*!
     * \brief Subtract the mean of constraint force from the background Eulerian
     * grid.
     */
    virtual void subtractMeanConstraintForce(Vec L, int f_data_idx, const double scale = 1.0);

    // \see CIBStrategy::setInterpolatedVelocityVector() method.
    /*!
     * \brief Prepare the CIBFEMethod class to get the interpolated fluid
     * velocity.
     */
    virtual void setInterpolatedVelocityVector(Vec V, const double data_time);

    // \see CIBStrategy::setInterpolatedVelocityVector() method.
    /*!
     * \brief Get interpolated velocity from the Eulerian grid.
     */
    virtual void getInterpolatedVelocity(Vec V, const double data_time, const double scale = 1.0);

    // \see CIBStrategy::computeMobilityRegularization method.
    /*!
     * \brief Compute regularization vector for the mobility problem.
     *
     */
    virtual void computeMobilityRegularization(Vec D, Vec L, const double scale = 1.0);

    // \see CIBStrategy::getNumberOfNodes method.
    /*!
     * \brief Get the number of nodes associated with a particular structure.
     */
    virtual unsigned int getNumberOfNodes(const unsigned int part) const;

    // \see CIBStrategy::computeNetRigidGeneralizedForce() methods.
    /*!
     * \brief Compute total force and torque on the rigid structure(s).
     */
    virtual void computeNetRigidGeneralizedForce(const unsigned int part, Vec L, RigidDOFVector& F);

    // \see CIBStrategy::copyVecToArray() method.
    /*!
     * \brief Copy PETSc Vec to raw array for specified structures.
     */
    void copyVecToArray(Vec b,
                        double* array,
                        const std::vector<unsigned>& struct_ids,
                        const int data_depth,
                        const int array_rank);

    // \see CIBStrategy::copyArrayToVec() method.
    /*!
     * \brief Copy raw array to PETSc Vec for specified structures.
     */
    void copyArrayToVec(Vec b,
                        double* array,
                        const std::vector<unsigned>& struct_ids,
                        const int data_depth,
                        const int array_rank);

    // \see CIBStrategy::setRigidBodyVelocity()
    /*!
     * \brief Set the rigid body velocity at the nodal points
     * contained in the Vec V.
     *
     */
    virtual void setRigidBodyVelocity(const unsigned int part, const RigidDOFVector& U, Vec V);

    // \}

    /*!
     * Initialize FE data.  This method must be called prior to calling
     * IBHierarchyIntegrator::initializePatchHierarchy().
     */
    void initializeFEData();

    /*!
     * \brief Register Eulerian variables with the parent IBHierarchyIntegrator.
     */
    virtual void registerEulerianVariables();

    /*!
     * \brief Register Eulerian refinement or coarsening algorithms with the parent
     * IBHierarchyIntegrator.
     */
    virtual void registerEulerianCommunicationAlgorithms();

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
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time);

    /*!
     * \brief Callbacks before INS is integrated.
     */
    typedef void (*preprocessSolveFluidEqn_callbackfcn)(const double, const double, const int, void*);

    /*!
     * \brief Register any preprocess fluid solve callback functions.
     */
    void registerPreProcessSolveFluidEquationsCallBackFcn(preprocessSolveFluidEqn_callbackfcn callback, void* ctx);

    /*!
     * \brief Calculate any body forces for INS solver over here.
     */
    virtual void preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num);

    /*!
     * Register a visIt data writer to output data files that
     * may be postprocessed with the visIt visualization tool.
     */
    void registerVisItDataWriter(SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * \brief Get the level number on which structures reside.
     */
    int getStructuresLevelNumber();

    /*!
     * \brief Get the patch hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > getPatchHierarchy();

    /*!
     *\brief Get the HierarchyMathOps object.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> getHierarchyMathOps();

    /*!
     * \brief Override the putToDatabase method of the base Serializable class.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    //////////////////////////////////////////////////////////////////////////////
protected:
    //////////////////////////////////////////////////////////////////////////////
private:
    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * \brief Get values from restart file.
     */
    void getFromRestart();

    /*!
     * \brief Compute the center of mass of the structure.
     */
    void computeCOMOfStructure(Eigen::Vector3d& center_of_mass, libMesh::EquationSystems* equation_systems);

    /*!
     * Number of nodes of rigid structures.
     */
    std::vector<unsigned int> d_num_nodes;

    /*!
     * FE data vectors.
     */
    std::vector<libMesh::System*> d_U_constrained_systems;
    std::vector<libMesh::PetscVector<double> *> d_U_constrained_current_vecs, d_U_constrained_half_vecs;
    std::vector<libMesh::PetscVector<double> *> d_F_current_vecs, d_F_new_vecs;

    /*!
     * Booleans to control spreading constraint force and interpolating
     * to Lagrangian velocities.
     */
    bool d_constraint_force_is_initialized, d_lag_velvec_is_initialized;

    /*!
     * Whether or not L2-projection is to be performed after velocity interpolation.
     */
    bool d_compute_L2_projection;

    /*!
     * Check if the initial center of mass has been initialized.
     */
    bool d_initial_com_initialized;

    /*!
     * PETSc wrappers for rigid body force.
     */
    std::vector<Vec> d_vL_new, d_vL_current;
    Vec d_mv_L_new, d_mv_L_current;

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
     * Control printing of S[lambda].
     */
    bool d_output_eul_lambda;

    /*!
     * Fluid density.
     */
    double d_rho;

}; // CIBFEMethod
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_CIBFEMethod
