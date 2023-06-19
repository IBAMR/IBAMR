// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
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

class LevelSetUtilities
{
public:
    /*!
     * \brief A lightweight class to hold the level set variable and the associated hierarchy
     * integrator (AdvDiffHierarchyIntegrator).
     */
    class LevelSetContainer : public virtual SAMRAI::tbox::DescribedClass
    {
    public:
        /*!
         * \brief Constructor of the class.
         */
        LevelSetContainer(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var)
            : d_adv_diff_integrator(adv_diff_integrator), d_ls_var(ls_var)
        {
            // return;
        }

        void setInterfaceHalfWidth(double ncells)
        {
            d_ncells = ncells;
        } // setInterfaceHalfWidth

        double getInterfaceHalfWidth()
        {
            return d_ncells;
        } // getInterfaceHalfWidth

        SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> getAdvDiffHierarchyIntegrator()
        {
            return d_adv_diff_integrator;
        } // getAdvDiffHierarchyIntegrator

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > getLevelSetVariable()
        {
            return d_ls_var;
        } // getLSVariable

    protected:
        SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_integrator;
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;
        double d_ncells = 1.0;
    };

    /*!
     * \brief A lightweight class to tag grid cells containing the level set variable for grid refinement
     */
    class TagLSRefinementCells : public LevelSetContainer
    {
    public:
        /*!
         * \brief Constructor of the class.
         */
        TagLSRefinementCells(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                             double tag_min_value = 0.0,
                             double tag_max_value = 0.0)
            : LevelSetContainer(adv_diff_integrator, ls_var),
              d_tag_min_value(tag_min_value),
              d_tag_max_value(tag_max_value)
        {
            return;
        } // TagLevelSetRefinementCells

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
        double d_tag_min_value = 0.0;
        double d_tag_max_value = 0.0;
    }; // TagLevelSetRefinementCells

    static void TagLSCells(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
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
    class LevelSetMassLossFixer : public LevelSetContainer, public SAMRAI::tbox::Serializable
    {
    public:
        /*!
         * \brief Constructor of the class.
         */
        LevelSetMassLossFixer(std::string object_name,
                              SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db =
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(nullptr),
                              bool register_for_restart = true);

        ~LevelSetMassLossFixer();

        void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

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
        } // getLagrangeMultiplier

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

    protected:
        std::string d_object_name;
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
     * \brief Compute the value of the Lagrange multiplier and use that to adjust the level set variable.
     */
    static void fixLevelSetMassLoss(double current_time,
                                    double new_time,
                                    bool skip_synchronize_new_state_data,
                                    int num_cycles,
                                    void* ctx);

    /*!
     * \return Integrals of the smooth Heaviside complement \f$ H(-phi) = 1- H(\phi)\f$ and smooth Heaviside function
     * \f$ H(\phi)\f$ over the entire domain.
     */
    static std::pair<double, double> computeIntegralHeavisideFcns(LevelSetContainer* lsc);

    /*!
     * \return Integral of inflow \f$ -\vec{u} \cdot \vec{n}\f$ at a physical boundary.
     *
     * \param location_idx of the boundary at which the integral is evaluated
     */
    static double computeNetInflowPhysicalBoundary(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                   int u_idx,
                                                   int bdry_loc_idx);

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
    static void setLSDataPatchHierarchy(int ls_idx,
                                        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                        const int integrator_step,
                                        const double current_time,
                                        const bool initial_time,
                                        const bool regrid_time,
                                        void* ctx);
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_LevelSetUtilities