// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_VelocityBcCoefs
#define included_VelocityBcCoefs

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/INSHierarchyIntegrator.h>

#include <ibtk/ibtk_utilities.h>

#include "BcData.h"

#include <RobinBcCoefStrategy.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class VelocityBcCoefs is an implementation of the strategy class
 * RobinBcCoefStrategy that is used to specify velocity boundary conditions that
 * are determined by a circulation model.
 */
class VelocityBcCoefs : public RobinBcCoefStrategyNd
{
public:
    /*!
     * \brief Constructor
     */
    VelocityBcCoefs(const INSHierarchyIntegrator* fluid_solver, const BcData, const int comp_idx);
    //    VelocityBcCoefs(const CirculationModel* circ_model, const int comp_idx);

    /*!
     * \brief Destructor.
     */
    virtual ~VelocityBcCoefs();

    /*!
     * \name Implementation of RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     */
    void setBcCoefs(Pointer<ArrayDataNd<double> >& acoef_data,
                    Pointer<ArrayDataNd<double> >& bcoef_data,
                    Pointer<ArrayDataNd<double> >& gcoef_data,
                    const Pointer<hier::VariableNd>& variable,
                    const PatchNd& patch,
                    const BoundaryBoxNd& bdry_box,
                    double fill_time) const;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     */
    IntVectorNd numberOfExtensionsFillable() const;

    //\}

    double parabolic_flow(double t, double y, double r) const;
    double time_ramp(double t) const;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VelocityBcCoefs(const VelocityBcCoefs& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VelocityBcCoefs& operator=(const VelocityBcCoefs& that);

    //    const CirculationModel* const d_circ_model;
    const INSHierarchyIntegrator* const d_fluid_solver;
    const int d_comp_idx;
    const BcData d_bc_data;
};

/////////////////////////////// INLINE ///////////////////////////////////////

// #include <VelocityBcCoefs.I>

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_VelocityBcCoefs
