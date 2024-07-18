// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_PhaseChangeDivUSourceFunction
#define included_IBAMR_PhaseChangeDivUSourceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/config.h>

#include "ibtk/CartGridFunction.h"

#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class PatchLevel;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

namespace IBAMR
{
class PhaseChangeHierarchyIntegrator;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class PhaseChangeDivUSourceFunction class set the RHS of the Div U equation computed from the
 * PhaseChangeHierarchyIntegrator class.
 */
class PhaseChangeDivUSourceFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    PhaseChangeDivUSourceFunction(const std::string& object_name,
                                  IBTK::SAMRAIPointer<IBAMR::PhaseChangeHierarchyIntegrator> pc_hier_integrator);

    /*!
     * \brief Empty destructor.
     */
    ~PhaseChangeDivUSourceFunction() = default;

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
                                 IBTK::SAMRAIPointer<SAMRAI::hier::VariableNd> var,
                                 IBTK::SAMRAIPointer<SAMRAI::hier::PatchHierarchyNd> hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        IBTK::SAMRAIPointer<SAMRAI::hier::VariableNd> var,
                        IBTK::SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                        const double data_time,
                        const bool initial_time = false,
                        IBTK::SAMRAIPointer<SAMRAI::hier::PatchLevelNd> patch_level =
                            IBTK::SAMRAIPointer<SAMRAI::hier::PatchLevelNd>(nullptr)) override;

private:
    /*!
     * Deleted default constructor.
     */
    PhaseChangeDivUSourceFunction() = delete;

    /*!
     * Deleted copy constructor.
     */
    PhaseChangeDivUSourceFunction(const PhaseChangeDivUSourceFunction& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    PhaseChangeDivUSourceFunction& operator=(const PhaseChangeDivUSourceFunction& that) = delete;

    /*!
     * Name of this object.
     */
    IBTK::SAMRAIPointer<IBAMR::PhaseChangeHierarchyIntegrator> d_pc_hier_integrator;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_PhaseChangeDivUSourceFunction
