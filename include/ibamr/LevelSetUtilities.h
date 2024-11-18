// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_LevelSetUtilities
#define included_IBAMR_LevelSetUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/LSInitStrategy.h"

#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include <limits>
#include <string>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
}
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
} // namespace pdat
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

namespace IBAMR
{
class AdvDiffHierarchyIntegrator;
} // namespace IBAMR

namespace IBTK
{
class HierarchyMathOps;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LevelSetUtilities provides implementation of some helper functions
 * to enforce mass/volume conservation of phases.
 */

namespace LevelSetUtilities
{
/*!
 * \brief A lightweight class to hold the level set variable and the associated hierarchy
 * integrator (AdvDiffHierarchyIntegrator).
 */
class LevelSetContainer : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor of the class.
     */
    LevelSetContainer(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                      double ncells = 1.0)
        : d_adv_diff_integrator(adv_diff_integrator), d_ncells(ncells)
    {
        d_ls_vars.emplace_back(ls_var);
        return;
    } // LevelSetContainer

    LevelSetContainer(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                      std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > ls_vars,
                      double ncells = 1.0)
        : d_adv_diff_integrator(adv_diff_integrator), d_ls_vars(std::move(ls_vars)), d_ncells(ncells)
    {
        // intentionally left blank
    } // LevelSetContainer

    /*!
     * @param ncells are the number of cells representing the half-width of the
     * interface. Depending upon the application, ncells is used to specify material
     * properties like density or viscosity (via some mixture model) or to compute volume
     * of the phase enclosed by the level set variable.
     */
    void setInterfaceHalfWidth(double ncells)
    {
        d_ncells = ncells;
    } // setInterfaceHalfWidth

    double getInterfaceHalfWidth() const
    {
        return d_ncells;
    } // getInterfaceHalfWidth

    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> getAdvDiffHierarchyIntegrator() const
    {
        return d_adv_diff_integrator;
    } // getAdvDiffHierarchyIntegrator

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > getLevelSetVariable(int idx = 0) const
    {
        return d_ls_vars[idx];
    } // getLSVariable

private:
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_integrator;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > d_ls_vars;
    double d_ncells = 1.0;
};

/*!
 * \brief A lightweight class to tag grid cells containing the level set variable for grid refinement
 */
class TagLSRefinementCells
{
public:
    /*!
     * \brief Constructor of the class.
     */
    TagLSRefinementCells(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                         double tag_min_value = 0.0,
                         double tag_max_value = 0.0)
        : d_ls_container(adv_diff_integrator, ls_var), d_tag_min_value(tag_min_value), d_tag_max_value(tag_max_value)
    {
        // intentionally left blank
    } // TagLevelSetRefinementCells

    const LevelSetContainer& getLevelSetContainer() const
    {
        return d_ls_container;
    } // getLevelSetContainer

    LevelSetContainer& getLevelSetContainer()
    {
        return d_ls_container;
    } // getLevelSetContainer

    void setTagMinValue(double tag_min)
    {
        d_tag_min_value = tag_min;
        return;
    } // setTagMinValue

    double getTagMinValue() const
    {
        return d_tag_min_value;
    } // getTagMinValue

    void setTagMaxValue(double tag_max)
    {
        d_tag_max_value = tag_max;
        return;
    } // setTagMaxValue

    double getTagMaxValue() const
    {
        return d_tag_max_value;
    } // getTagMaxValue

private:
    LevelSetContainer d_ls_container;
    double d_tag_min_value = 0.0;
    double d_tag_max_value = 0.0;
}; // TagLevelSetRefinementCells

/*!
 * \brief Preprocessing call back function to be hooked into IBAMR::HierarchyIntegrator class
 * to tag the cells for grid refinement based on the given tagging criteria.
 *
 * This static member should be registered with an appropriate hierarchy integrator
 * via registerApplyGradientDetectorCallback().
 *
 * \param ctx is the pointer to the TagLSRefinementCells class object.
 */
void tagLSCells(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                const int level_number,
                const double error_data_time,
                const int tag_index,
                const bool initial_time,
                const bool uses_richardson_extrapolation_too,
                void* ctx);

/*!
 * \brief A lightweight class that stores the current value of the Lagrange multiplier
 *  for the level set variable.
 */
class LevelSetMassLossFixer : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Constructor of the class.
     * @param input_db provides parameters such as enable_logging, correction_interval, max_its, rel_tol, and
     * half_width.
     */
    LevelSetMassLossFixer(std::string object_name,
                          SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                          std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > > ls_vars,
                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr,
                          bool register_for_restart = true);

    ~LevelSetMassLossFixer();

    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    const LevelSetContainer& getLevelSetContainer() const
    {
        return d_ls_container;
    } // getLevelSetContainer

    LevelSetContainer& getLevelSetContainer()
    {
        return d_ls_container;
    } // getLevelSetContainer

    // Getter and setter functions
    void setInitialVolume(double v0);

    double getInitialVolume() const
    {
        return d_vol_init;
    } // getInitialVolume

    void setTargetVolume(double v)
    {
        d_vol_target = v;
        return;
    } // setTargetVolume

    double getTargetVolume() const
    {
        return d_vol_target;
    } // setTargetVolume

    void setLagrangeMultiplier(double q)
    {
        d_q = q;
        return;
    } // setLagrangeMultiplier

    double getLagrangeMultiplier() const
    {
        return d_q;
    } // getLagrangeMultiplier

    void setTime(double time)
    {
        d_time = time;
        return;
    } // setTime

    double getTime() const
    {
        return d_time;
    } // getTime

    int getCorrectionInterval() const
    {
        return d_interval;
    } // getCorrectionInterval

    double getErrorRelTolerance() const
    {
        return d_rel_tol;
    } // getErrorRelTolerance

    double getMaxIterations() const
    {
        return d_max_its;
    } // getMaxIterations

    bool enableLogging() const
    {
        return d_enable_logging;
    } // enableLogging

private:
    std::string d_object_name;
    LevelSetContainer d_ls_container;
    bool d_registered_for_restart, d_enable_logging = false;

    /*
     * Initial volume of the phase that is enforced by the Lagrange multiplier
     */
    double d_vol_init = std::numeric_limits<double>::quiet_NaN();

    /*
     * Target volume of the phase that is enforced by the Lagrange multiplier
     */
    double d_vol_target = std::numeric_limits<double>::quiet_NaN();

    /*
     * The Lagrange multiplier which enforces volume conservation.
     */
    double d_q = std::numeric_limits<double>::quiet_NaN();

    double d_time = 0.0;

    int d_interval = 1, d_max_its = 4;

    double d_rel_tol = 1e-12;

    void getFromRestart();

    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);
};

/*!
 * \brief Compute the value of the Lagrange multiplier \f$ q \f$ and use that to adjust the flow level set variable
 * \f$ \tilde{\phi} \f$ to satisfy the constraint: \f$ f(q) = \int_{\Omega} H(\tilde{\phi} + q) \text{d}\Omega - V^0
 * = 0 \f$. where, \f$ \tilde{\phi} \f$ is the level set field obtained from the reinitialization procedure, \f$ V^0
 * \f$ is the volume of the fluid at \f$ t = 0 \f$ s which needs to be conserved. Here, \f$ q \f$ is computed using
 * the Newton's method till required tolerance. In practice the level set mass loss would be fixed during the
 * postprocess integrate hierarchy stage. Hence the application time would be the new time and the variable context
 * would be the new context.
 *
 * \param ctx is the pointer to the LevelSetMassLossFixer class object.
 */
void fixMassLoss2PhaseFlows(double current_time,
                            double new_time,
                            bool skip_synchronize_new_state_data,
                            int num_cycles,
                            void* ctx);

/*!
 * \brief Compute the value of the Lagrange multiplier \f$ q \f$ and use that to adjust the flow level set \f$
 * \tilde{\phi \f$ when there is a solid phase in the domain with level set \f$ \psi < 0 \f$. Satisfies the
 * constraint: \f$ f(q) = \int_{\Omega} H(\tilde{\phi} + q) H(\psi) \text{d}\Omega - V^0 = 0 \f$. where, \f$
 * \tilde{\phi} \f$ is the level set field obtained from the reinitialization procedure, \f$ V^0 \f$ is the volume
 * of the fluid at \f$ t = 0 \f$ s which needs to be conserved. Here, \f$ q \f$ is computed using the Newton's
 * method till required tolerance. In practice the level set mass loss would be fixed during the postprocess
 * integrate hierarchy stage. Hence the application time would be the new time and the variable context would be the
 * new context.
 *
 * \param ctx is the pointer to the LevelSetMassLossFixer class object.
 */
void fixMassLoss3PhaseFlows(double current_time,
                            double new_time,
                            bool skip_synchronize_new_state_data,
                            int num_cycles,
                            void* ctx);

/*!
 * \return Integrals of \f$ 1- H(\phi)\f$, \f$ H(\phi)\f$, and \f$ \delta(\phi) \f$  over
 * the entire domain.
 *
 * Here, \f$ H(\phi) \f$ is the Heaviside function demarcating liquid and gas domains.
 *
 * \f$ \phi \f$ is taken to be positive in the liquid domain and negative in the gas domain.
 *
 * Physically, these three integrals represent volume of the gas region, liquid region, and
 * surface area of the interface, respectively.
 */
std::vector<double> computeHeavisideIntegrals2PhaseFlows(const LevelSetContainer& lsc);

/*!
 * \return Integrals of \f$ [1- H(\phi)] H(\Psi) \f$, \f$ H(\phi) H(\Psi)\f$, \f$ 1 - H(\Psi)\f$,
 * and \f$ \delta(\phi)H(\Psi)\f$ over the entire domain.
 *
 * Here, \f$ H(\phi) \f$ is the Heaviside function demarcating liquid and gas domains,  and \f$ H(\Psi) \f$ is
 * the Heaviside function demarcating solid and fluid (fluid = liquid and gas) domains.
 *
 * Physically, these four integrals represents volume of the gas region, liquid region, solid region, and
 * surface area of the fluid interface, respectively.
 *
 * \f$ \phi \f$ is taken to be positive in the liquid domain and negative in the gas domain.
 *
 * \f$ \Psi \f$ is taken to be positive outside the solid and negative inside the solid.
 */
std::vector<double> computeHeavisideIntegrals3PhaseFlows(const LevelSetContainer& lsc);

/*!
 * \brief Class SetLSProperties is a utility class which sets (or resets after reinitialization)
 * level set values on the patch hierarchy.
 */
class SetLSProperties
{
public:
    /*!
     * The only constructor of this class.
     */
    SetLSProperties(const std::string& object_name, SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> ls_ops)
        : d_object_name(object_name), d_ls_ops(ls_ops)
    {
        // intentionally left blank
        return;
    } // SetLSProperties

    /*!
     * Destructor for this class.
     */
    ~SetLSProperties() = default;

    /*!
     * Set the level set value on the patch hierarchy
     */
    void setLSData(int ls_idx,
                   SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                   const int integrator_step,
                   const double current_time,
                   const bool initial_time,
                   const bool regrid_time);

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetLSProperties() = delete;

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetLSProperties& operator=(const SetLSProperties& that) = delete;

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetLSProperties(const SetLSProperties& from) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    SAMRAI::tbox::Pointer<IBAMR::LSInitStrategy> d_ls_ops;

}; // SetLSProperties

/*!
 * \brief A function that sets or resets the level set data on the patch hierarchy
 */
void setLSDataPatchHierarchy(int ls_idx,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int integrator_step,
                             const double current_time,
                             const bool initial_time,
                             const bool regrid_time,
                             void* ctx);
} // namespace LevelSetUtilities

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_LevelSetUtilities
