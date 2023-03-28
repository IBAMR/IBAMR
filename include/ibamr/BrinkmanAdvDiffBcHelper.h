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

#ifndef included_BrinkmanAdvDiffBcHelper
#define included_BrinkmanAdvDiffBcHelper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/config.h"

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/ibtk_utilities.h"

#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
class HierarchyMathOps;
}
/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief BrinkmanAdvDiffBcHelper is an abstract class that provides an interface
 * to implement Brinkman penalization body force in the advection-diffusion equation
 * in order to enforce Dirichlet, Neumann and Robin boundary conditions on surfaces of rigid
 * immersed bodies. A single instance of this class is meant to handle all of the Brinkman
 * penalization zones for multiple transported quantities with various boundary conditions.

 * BrinkmanAdvDiffBcHelper provides an implementation of a volume penalized
 * body force and linear operator modifications required to impose Dirichlet, Neumann and Robin
 * boundary conditions to scalar quantities maintained by BrinkmanSemiImplicitAdvDiffHierarchyIntegrator.
 *
 * Boundary conditions can be applied to multiple interfaces, which are demarcated
 * using level set variables. This class assumes that the penalized region coincides with negative
 * values of the level set. The sign convention of the level set variable is specified by the user.
 *
 * Reference
 * Sakurai, T., Yoshimatsu, K., Okamoto N. and Schneider K.,<A
 HREF="https://www.sciencedirect.com/science/article/pii/S0021999119302414">
 * Volume penalization for inhomogeneous Neumann boundary conditions modeling scalar flux in complicated geometry</A>
 */
class BrinkmanAdvDiffBcHelper : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor of the class.
     */
    BrinkmanAdvDiffBcHelper(std::string object_name,
                            SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver);

    /*!
     * \brief Destructor of the class.
     */
    ~BrinkmanAdvDiffBcHelper() = default;

    /*!
     * \brief Set the time interval in which Brinkman forcing is computed.
     */
    void setTimeInterval(double current_time, double new_time);

    /*!
     * \brief Preprocess routine before computing Brinkman penalization terms.
     *
     */
    void preprocessBrinkmanAdvDiffBcHelper(double current_time, double new_time, int num_cycles);

    /*!
     * \brief Postprocess routine after computing Brinkman penalization terms.
     *
     */
    void postprocessBrinkmanAdvDiffBcHelper(double current_time, double new_time, int num_cycles);

    /*!
     * \brief Set Brinkman penalization penalty factor for all level sets.
     */
    void setPenaltyCoefficient(double eta_penalty_coef);

    /*!
     * \brief Set the number of interface cells for all level sets.
     */
    void setNumInterfaceCells(double num_interface_cells);

    /*!
     * \brief Get the name of the object.
     */
    const std::string& getName() const
    {
        return d_object_name;
    } // getName

    /*!
     * \brief Get the current time interval \f$ [t^{n+1}, t^n] \f$ in which Brinkman
     * velocity is computed.
     */
    std::pair<double, double> getCurrentTimeInterval() const
    {
        return std::make_pair(d_new_time, d_current_time);
    } // getCurrentTimeInterval

    /*!
     * \brief Function to determine if a transported quantity has Brinkman
     * boundary conditions associated with it.
     */
    bool hasBrinkmanBoundaryCondition(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var) const
    {
        return d_Q_bc.find(Q_var) != d_Q_bc.end();
    } // hasBrinkmanBoundaryCondition

    /*!
     * \brief Function specifying the optional forcing function for inhomogeneous
     * boundary conditions \f$ \zeta q + \kappa n \dot \nabla q = \g \f$
     *
     * The user must set the patch data B_idx such that \f$ n \dot B = g \f$.
     * The applied forcing term is then computed internally as \f$ \nabla \dot (\chi B) - \chi \nabla \dot B \f$.
     * Note that B_idx contains side-centered patch data.
     */
    using BrinkmanInhomogeneousBCsFcnPtr =
        void (*)(int B_idx,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                 SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                 double time,
                 void* ctx);

    /*!
     * \brief Register a transported quantity with this object, along with the solid level set
     * variable for which to apply a homogeneous boundary condition.
     *
     * \note This function can only be used to register homogeneous Dirichlet, Neumann and Robin
     * BCs.
     */
    void registerHomogeneousBC(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                               SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                               std::string bc_type,
                               std::string indicator_func_type = "SMOOTH",
                               double num_interface_cells = 2.0,
                               double eta_penalty_coef = 1.0e-8);

    /*!
     * \brief Register a transported quantity with this object, along with the solid level set
     * variable on which to apply inhomogeneous boundary
     * conditions.
     *
     * \note Inhomogeneous BCs are treated uniquely within this class and require
     * additional user callback inputs.
     */
    void registerInhomogeneousBC(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                                 std::string bc_type,
                                 BrinkmanInhomogeneousBCsFcnPtr callback = nullptr,
                                 void* ctx = nullptr,
                                 std::string indicator_func_type = "SMOOTH",
                                 double num_interface_cells = 2.0,
                                 double eta_penalty_coef = 1.0e-8,
                                 double bc_val = std::numeric_limits<double>::signaling_NaN());

    /*!
     * \brief Function to compute the cell-centered coefficient to the damping linear operator
     * and RHS of the advection-diffusion equation for a specified transported quantity Q_var
     *
     * \note It is assumed that the physical damping coefficient \f$\lambda\f$ is zero.
     *
     * The functional form of the Brinkman damping coefficient is
     * \f$ C = \sum_{Dirichlet} \chi_i/\eta  + \sum_{Robin} \nabla \dot (\chi_i n_i) - \chi_i \nabla \dot n_i\f$
     * where the sum is taken over all level sets with Dirichlet and Robin BCs. Here \f$\chi_i = 1-H_i\f$. Note that
     * Neumann BCs do not contribute anything to this term.
     *
     */
    void computeDampingCoefficient(int C_idx, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * \brief Function to compute the side-centered coefficient to the diffusion linear operator
     * and RHS of the advection-diffusion equation for a specified transported quantity Q_var with diffusion
     * coefficient kappa.
     *
     * \note This function is able to handle both constant and variable kappa.
     *
     * The functional form of the Brinkman diffusion coefficient is
     * \f$ D = \kappa + \sum_{Neumann} (-\chi_i + \eta \chi_i) + \sum_{Robin} (-\chi_i + \eta \chi_i)\f$
     * where \f$\chi_i = 1-H_i\f$ and the sum is taken over
     * all level sets with Neumann and Robin BCs . Note that Dirichlet BCs do not contribute anything to this term.
     */
    void computeDiffusionCoefficient(int D_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                     int kappa_idx,
                                     double kappa);

    /*!
     * \brief Function to compute the Brinkman forcing contribution to the RHS of the advection-diffusion solver
     * for a specified transported quantity Q_var.
     *
     * For Inhomogeneous Dirichlet BCs, \f$ F_i = \chi_i/\eta Q_{bc}\f$ where \f$\chi_i = 1-H_i\f$.
     * For inhomogeneous Neumann and Robin BCs, \f$ F_i = \nabla \dot (\chi_i B) - \chi_i \nabla \dot B_i \f$,
     * with a user defined \f$B_i\f$.
     * The overall functional form of the Brinkman body force is \f$F = \sum_{i} F_i\f$.
     */
    void computeForcing(int F_idx, SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var);

    /*!
     * \brief Function to mask the additional forcing terms on the RHS of the advection-diffusion solver
     * e.g. \f$ u \dot \grad Q\f$ and body forces.
     *
     * The functional form of the Brinkman masking term is
     * \f$ N = (1-\sum_{Neumann and Robin} \chi_i) N\f$ where \f$\chi_i = 1-H_i\f$ and the sum
     * is taken over all level sets with Neumann and Robin BCs. Note that Dirichlet BCs do
     * not mask this term presently.
     */
    void maskForcingTerm(int N_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                         const bool mask_smeared_region = false);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BrinkmanAdvDiffBcHelper(const BrinkmanAdvDiffBcHelper& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BrinkmanAdvDiffBcHelper& operator=(const BrinkmanAdvDiffBcHelper& that) = delete;

    /*!
     * Struct to maintain the properties a particular boundary condition
     */
    struct BCProperties
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var;
        AdvDiffBrinkmanPenalizationBcType bc_type;
        double bc_val;
        double num_interface_cells;
        double eta;
        IndicatorFunctionType indicator_func_type;
        BrinkmanInhomogeneousBCsFcnPtr callback;
        void* ctx;
    };

    /*!
     * Object name.
     */
    std::string d_object_name;

    /*!
     * Pointer to the adv-diff solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Time interval
     */
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN();

    /*!
     * Map for storing transported variables and various registered BCs
     */
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >, std::vector<BCProperties> > d_Q_bc;

    /*!
     * Patch data required for computing additional forcing for inhomogeneous Neumann and Robin BCs
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_B_var;
    int d_B_scratch_idx = IBTK::invalid_index, d_B_chi_scratch_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_div_var;
    int d_div_B_scratch_idx = IBTK::invalid_index, d_div_B_chi_scratch_idx = IBTK::invalid_index,
        d_chi_scratch_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_n_var;
    int d_n_scratch_idx = IBTK::invalid_index, d_n_chi_scratch_idx = IBTK::invalid_index;
    int d_div_n_scratch_idx = IBTK::invalid_index, d_div_n_chi_scratch_idx = IBTK::invalid_index;
    int d_variable_g_scratch_idx = IBTK::invalid_index;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_BrinkmanAdvDiffBcHelper
