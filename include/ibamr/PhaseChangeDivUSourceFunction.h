// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2018 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_PhaseChangeDivUSourceFunction
#define included_PhaseChangeDivUSourceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/PhaseChangeHierarchyIntegrator.h>

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class PhaseChangeDivUSourceFunction class set the RHS of the Div U equation computed from the
 * PhaseChangeHierarchyIntegrator class.
 */
class PhaseChangeDivUSourceFunction : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    PhaseChangeDivUSourceFunction(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<IBAMR::PhaseChangeHierarchyIntegrator> pc_hier_integrator,
                                  SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops);

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
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL)) override;

    //\}

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
    std::string d_object_name;

    /*!
     * Name of this object.
     */
    SAMRAI::tbox::Pointer<IBAMR::PhaseChangeHierarchyIntegrator> d_pc_hier_integrator;

    /*!
     * Pointer to HierarchyMathOps.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PhaseChangeDivUSourceFunction
