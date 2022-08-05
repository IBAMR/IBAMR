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

#ifndef included_LevelSetInitialCondition
#define included_LevelSetInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibtk/muParserCartGridFunction.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class LevelSetInitialCondition provides an initial condition for
 * the level set function.
 */
class LevelSetInitialCondition : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LevelSetInitialCondition(const std::string& object_name, const double initial_gas_pcm_interface_position);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetInitialCondition() = default;

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
    LevelSetInitialCondition() = delete;

    /*!
     * Deleted copy constructor.
     */
    LevelSetInitialCondition(const LevelSetInitialCondition& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    LevelSetInitialCondition& operator=(const LevelSetInitialCondition& that) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Initial position of a gas-pcm interface.
     */
    const double d_initial_gas_pcm_interface_position;
};
//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LevelSetInitialCondition
