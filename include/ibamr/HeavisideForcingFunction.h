// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_HeavisideForcingFunction
#define included_HeavisideForcingFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include "ibtk/samrai_compatibility_names.h"
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>

#include "SAMRAICellVariable.h"
#include "SAMRAIFaceVariable.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIVariable.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBAMR
{
/*!
 * \brief Class HeavisideForcingFunction computes \f$ H \nabla \cdot \vec{u} \f$. This is required to solve
 * the advection equation for Heaviside in conservative form.
 */
class HeavisideForcingFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    HeavisideForcingFunction(const std::string& object_name,
                             SAMRAIPointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                             SAMRAIPointer<SAMRAICellVariable<double>> H_var,
                             SAMRAIPointer<SAMRAIFaceVariable<double>> U_adv_var);

    /*!
     * \brief Empty destructor.
     */
    ~HeavisideForcingFunction() = default;

    /*!
     * \brief Indicates whether the concrete TemperatureInitialCondition object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy using the virtual function
     * setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAIPointer<SAMRAIVariable> var,
                                 SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAIPointer<SAMRAIVariable> var,
                        SAMRAIPointer<SAMRAIPatch> patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAIPointer<SAMRAIPatchLevel> patch_level = nullptr) override;

    //\}

private:
    /*!
     * Deleted default constructor.
     */
    HeavisideForcingFunction() = delete;

    /*!
     * Deleted copy constructor.
     */
    HeavisideForcingFunction(const HeavisideForcingFunction& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    HeavisideForcingFunction& operator=(const HeavisideForcingFunction& that) = delete;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAIPointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Liquid fraction variable.
     */
    SAMRAIPointer<SAMRAICellVariable<double>> d_H_var;

    /*!
     * Advection velocity variable.
     */
    SAMRAIPointer<SAMRAIFaceVariable<double>> d_U_adv_var;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_HeavisideForcingFunction
