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

#ifndef included_LiquidFractionInitialCondition
#define included_LiquidFractionInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibtk/muParserCartGridFunction.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class LiquidFractionInitialCondition provides an initial condition for
 * the liquid fraction. In this example, the liquid and solid occupies the bottom and middle regions,
 * respectively. The gas region is placed on top of the solid and we set the liquid fraction to be zero
 * in the gas.
 */
class LiquidFractionInitialCondition : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LiquidFractionInitialCondition(const std::string& object_name,
                                   const double initial_liquid_solid_interface_position);

    /*!
     * \brief Empty destructor.
     */
    ~LiquidFractionInitialCondition() = default;

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
    LiquidFractionInitialCondition() = delete;

    /*!
     * Deleted copy constructor.
     */
    LiquidFractionInitialCondition(const LiquidFractionInitialCondition& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    LiquidFractionInitialCondition& operator=(const LiquidFractionInitialCondition& that) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Initial position of liquid-solid interface.
     */
    const double d_initial_liquid_solid_interface_position;
};
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LiquidFractionInitialCondition
