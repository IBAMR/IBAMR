// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_LevelSetGasInitialCondition
#define included_LevelSetGasInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class LevelSetGasInitialCondition provides an initial condition for
 * the level set function.
 */
class LevelSetGasInitialCondition : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LevelSetGasInitialCondition(const std::string& object_name,
                                const double greater_x_column,
                                const double less_z_column);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetGasInitialCondition();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete LevelSetGasInitialCondition object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL));

    //\}

private:
    LevelSetGasInitialCondition();

    LevelSetGasInitialCondition(const LevelSetGasInitialCondition& from);

    LevelSetGasInitialCondition& operator=(const LevelSetGasInitialCondition& that);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Initial level set information.
     */
    double d_greater_x_column, d_less_z_column;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LevelSetGasInitialCondition
