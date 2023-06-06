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

#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
} // namespace pdat
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
     * \brief A lightweight structure to hold the level set variable and the AdvDiffHierarchyIntegrator
     */
    struct LevelSetContainer
    {
        /*!
         * \brief Constructor of the class.
         */
        LevelSetContainer(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                          int ncells = 1)
            : d_adv_diff_integrator(adv_diff_integrator), d_ls_var(ls_var), d_ncells(ncells)
        {
            // return;
        }

        SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_integrator;
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;
        int d_ncells;
        double d_time;
    };

    /*!
     * \brief A lightweight structure that stores the current value of the Lagrange multiplier
     *  for the level set variable.
     */
    struct LevelSetMassLossFixer : public LevelSetContainer
    {
        /*!
         * \brief Constructor of the class.
         */
        LevelSetMassLossFixer(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                              int ncells = 1)
            : LevelSetContainer(adv_diff_integrator, ls_var, ncells)
        {
            return;
        } // LevelSetMassLossFixer

        double d_vol_init, d_q;
    };

    /*!
     * \brief Compute the value of the Lagrange multiplier and use that to adjust the level set variable.
     */
    static void fixLevelSetMassLoss(double current_time, double new_time, int cycle_num, void* ctx);

    /*!
     * \return Integral of the smooth Heaviside function \f$ H(\phi)\f$ and its complement \f$ H(-phi) = 1- H(\phi)\f$
     * over the entire domain.
     */
    static std::pair<double, double>
    computeIntegralHeavisideFcns(double current_time, double new_time, int cycle_num, void* ctx);

    /*!
     * \brief Class SetLSProperties is a utility class which sets level set values on the patch hierarchy
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
    static void setLSDataHierarchy(int ls_idx,
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