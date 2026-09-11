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
                                double initial_liquid_solid_interface_position,
                                double initial_liquid_temperature,
                                double initial_solid_temperature);

    /*!
     * \brief Empty destructor.
     */
    ~TemperatureInitialCondition() override = default;

    /*!
     * \brief Indicates whether the concrete TemperatureInitialCondition object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        Pointer<Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        Pointer<PatchLevel<NDIM>> patch_level = Pointer<PatchLevel<NDIM>>(nullptr)) override;

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
    double d_initial_liquid_solid_interface_position;

    /*!
     * Initial temperature of the liquid.
     */
    double d_initial_liquid_temperature;

    /*!
     * Initial temperature of the solid.
     */
    double d_initial_solid_temperature;
};
//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_TemperatureInitialCondition
