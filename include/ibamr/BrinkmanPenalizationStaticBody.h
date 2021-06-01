// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_BrinkmanPenalizationStaticBody
#define included_BrinkmanPenalizationStaticBody

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/BrinkmanPenalizationStrategy.h"
#include "ibamr/IBHydrodynamicSurfaceForceEvaluator.h"

#include "ibtk/ibtk_utilities.h"

#include "tbox/Pointer.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "Eigen/Core"
#include "Eigen/Geometry"
IBTK_ENABLE_EXTRA_WARNINGS

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
class INSVCStaggeredHierarchyIntegrator;
class AdvDiffHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief BrinkmanPenalizationStaticBody provides an implementation of Brinkman penalization
 * body force for a rigid body motion in the momentum equation.
 *
 * The penalization force is taken to be \f$ \frac{\chi}{\kappa}(\bm{u}_b
 * - \bm{u}^{n+1}) \f$. The class computes the coefficient \f$
 * \frac{\chi}{\kappa}$ of the fluid velocity \f$ \bm{u}^{n+1} \f$ for the
 * variable-coefficient INS solvers of INSVCStaggeredHierarchyIntegrator. This
 * is done in the BrinkmanPenalizationStaticBody::demarcateBrinkmanZone
 * method. Here \f$ \chi \f$ is the body indicator function and \f$\kappa \sim
 * \Delta t/ \rho \ll 1 \f$ is the vanishing permeability of the body. The
 * rigid body velocity \f$\bm{u}_b\f$ is computed through Newton's law of
 * motion by netting the hydrodynamic and external forces and torques on the
 * body in the BrinkmanPenalizationStaticBody::computeBrinkmanVelocity
 * method. A simple forward-Euler scheme is employed for the second-law of
 * motion.
 *
 * For further information on applications of this class see
 * https://arxiv.org/abs/1904.04078.
 */
class BrinkmanPenalizationStaticBody : public BrinkmanPenalizationStrategy
{
public:
    /*!
     * Since this class has Eigen object members, which have special alignment
     * requirements, we must explicitly override operator new to get the
     * correct aligment for the object as a whole.
     */
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /*
     * \brief Constructor of the class.
     */
    BrinkmanPenalizationStaticBody(std::string object_name,
                                   SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                                   SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                   SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                                   SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> fluid_solver,
                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr,
                                   bool register_for_restart = true);

    /*
     * \brief Destructor of the class.
     */
    ~BrinkmanPenalizationStaticBody() = default;

    /*!
     * \brief Preprocess routine before computing Brinkman penalization terms.
     *
     */
    void preprocessComputeBrinkmanPenalization(double current_time, double new_time, int num_cycles) override;

    /*!
     * \brief Compute the desired rigid body velocity in the Brinkman penalized zone.
     */
    void computeBrinkmanVelocity(int u_idx, double time, int cycle_num) override;

    /*!
     * \brief Demarcate the Brinkman zone with Brinkman penalty term.
     */
    void demarcateBrinkmanZone(int u_idx, double time, int cycle_num) override;

    /*!
     * \brief Postprocess routine after computing Brinkman penalization terms.
     */
    void postprocessComputeBrinkmanPenalization(double current_time, double new_time, int num_cycles) override;

    /*!
     * \brief Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    // Mass and inertial of the body.
    double d_mass = 0.0;
    Eigen::Matrix3d d_inertia_tensor_initial = Eigen::Matrix3d::Zero();

    // Pointers to solvers.
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_fluid_solver;

    // Level set variable defining the solid.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_solid_var;

    // Liquid fraction variable.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;

    // Hydrodynamic force evaluator.
    SAMRAI::tbox::Pointer<IBAMR::IBHydrodynamicSurfaceForceEvaluator> d_hydro_force_eval;

    // Contour level
    double d_contour_level = 0.0;

    // Number of interface cells to compute the Heaviside function
    double d_num_interface_cells = 2.0;

    // Forces and torques on the body.
    Eigen::Vector3d d_hydro_force_pressure, d_hydro_force_viscous, d_hydro_torque_pressure, d_hydro_torque_viscous,
        d_ext_force, d_ext_torque;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BrinkmanPenalizationStaticBody(const BrinkmanPenalizationStaticBody& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BrinkmanPenalizationStaticBody& operator=(const BrinkmanPenalizationStaticBody& that) = delete;

    /*!
     * \brief Get options from input database.
     */
    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Constant added to avoid zero division.
     */
    const double d_ed = 1e-3;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_BrinkmanPenalizationStaticBody
