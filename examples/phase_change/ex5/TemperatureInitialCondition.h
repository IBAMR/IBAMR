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

#ifndef included_TemperatureInitialCondition
#define included_TemperatureInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibtk/muParserCartGridFunction.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class TemperatureInitialCondition provides an initial condition for
 * the temperature. In this example, the liquid and solid occupies the bottom and middle regions,
 * respectively. The gas region is placed on top of the solid and we set the solid temperature in
 * the gas as an initial condition.
 */
class TemperatureInitialCondition : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    TemperatureInitialCondition(const std::string& object_name,
                                const double initial_liquid_solid_interface_position,
                                const double initial_liquid_temperature,
                                const double initial_solid_temperature);

    /*!
     * \brief Empty destructor.
     */
    ~TemperatureInitialCondition() = default;

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
    TemperatureInitialCondition() = delete;

    /*!
     * Deleted copy constructor.
     */
    TemperatureInitialCondition(const TemperatureInitialCondition& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    TemperatureInitialCondition& operator=(const TemperatureInitialCondition& that) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Initial position of liquid-solid interface.
     */
    const double d_initial_liquid_solid_interface_position;

    /*!
     * Initial temperature of the liquid.
     */
    const double d_initial_liquid_temperature;

    /*!
     * Initial temperature of the solid.
     */
    const double d_initial_solid_temperature;
};
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_TemperatureInitialCondition
