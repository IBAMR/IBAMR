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

#ifndef included_CarmanKozenyDragForce
#define included_CarmanKozenyDragForce

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/BrinkmanPenalizationStrategy.h"

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
 * \brief CarmanKozenyDragForce provides an implementation of Carman-Kozeny drag force
 * to impose zero velocity inside the solid \f$ \bm{u}=\bm{u}_b = 0\f$.
 *
 * The penalization force is taken to be \f$ A_d(\bm{u}_b - \bm{u}^{n+1}) \f$. The class
 * computes the coefficient \f$A_d \f$ of the fluid velocity \f$ \bm{u}^{n+1} \f$ for the
 * variable-coefficient INS solvers of INSVCStaggeredHierarchyIntegrator. This is done in
 * the CarmanKozenyDragForce::demarcateBrinkmanZone method. Here \f$ A_d =
 * C_d\frac{\alpha_S}{(1-alpha_S)^3+e_d}\f$, \f$ \alpha_S \f$ is the volume fraction of
 * the solid, \f$ C_d\f$ and  \f$ e_d\f$ are the model parameters. The rigid body velocity
 * \f$\bm{u}_b\f$ is taken to be zero. The penalty
 * parameter C_d is taken to be \f$C_d = ( \rho / \Delta t + \mu / h^2)\f$. The user can choose
 * the density or the inertia scale or both for \f$C_d\f$ through input.
 */
class CarmanKozenyDragForce : public BrinkmanPenalizationStrategy
{
public:
    /*
     * \brief Constructor of the class.
     */
    CarmanKozenyDragForce(std::string object_name,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > H_var,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                          SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                          SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> fluid_solver,
                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr,
                          bool register_for_restart = true);

    /*
     * \brief Destructor of the class.
     */
    ~CarmanKozenyDragForce() = default;

    /*!
     * \brief Preprocess routine before computing Carman-Kozeny term.
     *
     */
    void preprocessComputeBrinkmanPenalization(double current_time, double new_time, int num_cycles) override;

    /*!
     * \brief Compute the desired rigid body velocity in the Brinkman penalized (solid) zone. The present
+    * implementation sets rigid body velocity to be zero.
     */
    void computeBrinkmanVelocity(int u_idx, double time, int cycle_num) override;

    /*!
     * \brief Demarcate the Brinkman zone (solid) with Carman-Kozeny term term.
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
    /*!
     * \brief Pointers to solvers.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_fluid_solver;

    /*!
     * \brief Heaviside variable defining the gas-pcm interface.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_H_var;

    /*!
     * \brief Liquid fraction variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CarmanKozenyDragForce(const CarmanKozenyDragForce& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CarmanKozenyDragForce& operator=(const CarmanKozenyDragForce& that) = delete;

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
     * Constant needed to avoid division by zero. This constant also controls the strength
     * of the penalty factor.
     */
    double d_ed = 1e-3;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_CarmanKozenyDragForce
