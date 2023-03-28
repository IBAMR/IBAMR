// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_FeedbackForcer
#define included_FeedbackForcer

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/INSHierarchyIntegrator.h>

#include "BcData.h"

#include <Patch.h>
#include <Variable.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class FeedbackForcer is an implementation of the strategy class
 * CartGridFunction that is used to specify velocity boundary conditions via a
 * feedback forcing (i.e., penalty) method.
 */
class FeedbackForcer : public CartGridFunction
{
public:
    /*!
     * \brief Constructor
     */
    FeedbackForcer(const INSHierarchyIntegrator* fluid_solver,
                   Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                   const BcData& bc_data);

    /*!
     * \brief Destructor.
     */
    virtual ~FeedbackForcer();

    /*!
     * \name Implementation of CartGridFunction interface.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Set data on the specified patch interior.
     */
    void setDataOnPatch(int data_idx,
                        Pointer<hier::Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL));

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FeedbackForcer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FeedbackForcer(const FeedbackForcer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FeedbackForcer& operator=(const FeedbackForcer& that);

    const INSHierarchyIntegrator* const d_fluid_solver;
    Pointer<PatchHierarchy<NDIM> > d_patch_hierarchy;
    const BcData d_bc_data;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <FeedbackForcer.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FeedbackForcer
