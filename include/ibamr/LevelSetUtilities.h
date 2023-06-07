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

#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include <limits>
#include <string>

namespace SAMRAI
{
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
class LSInitStrategy;
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
         *
         * @param ncells The number of cells across either side of the interface to smear the Heaviside function.
         */
        LevelSetContainer(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                          double ncells = 1)
            : d_adv_diff_integrator(adv_diff_integrator), d_ls_var(ls_var), d_ncells(ncells)
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
        double d_ncells;
    };

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
                              double ncells = 1.0,
                              bool register_for_restart = true);

        ~LevelSetMassLossFixer();

        void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

        void setInitialVolume(double v0);

        double getInitialVolume() const
        {
            return d_vol_init;
        } // getInitialVolume

        void setLagrangeMultiplier(double q)
        {
            d_q = q;
        } // getLagrangeMultiplier

        double getLagrangeMultiplier() const
        {
            return d_q;
        } // getLagrangeMultiplier

        void setTime(double time)
        {
            d_time = time;
        } // setTime

        double getTime() const
        {
            return d_time;
        } // getTime

    protected:
        std::string d_object_name;
        bool d_registered_for_restart;

        /*
         * Initial volume of the phase that is enforced by the Lagrange multiplier
         */
        double d_vol_init = std::numeric_limits<double>::quiet_NaN();

        /*
         * The Lagrange multiplier which enforces volume conservation.
         */
        double d_q = std::numeric_limits<double>::quiet_NaN();

        double d_time = 0.0;

        void getFromRestart();
    };

    /*!
     * \brief Compute the value of the Lagrange multiplier and use that to adjust the level set variable.
     */
    static void fixLevelSetMassLoss(double current_time, double new_time, int cycle_num, void* ctx);

    /*!
     * \return Integrals of the smooth Heaviside function \f$ H(\phi)\f$ and its complement \f$ H(-phi) = 1- H(\phi)\f$
     * over the entire domain.
     */
    static std::pair<double, double>
    computeIntegralHeavisideFcns(double current_time, double new_time, int cycle_num, void* ctx);

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