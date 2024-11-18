// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_hagen_poiseuille_FeedbackForcer
#define included_hagen_poiseuille_FeedbackForcer

#include <ibamr/INSHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class hagen_poiseuille_FeedbackForcer is an implementation of the strategy class
 * CartGridFunction that is used to specify velocity boundary conditions via a
 * feedback forcing (i.e., penalty) method.
 */
class hagen_poiseuille_FeedbackForcer : public CartGridFunction
{
public:
    /*!
     * \brief Constructor
     */
    hagen_poiseuille_FeedbackForcer(const double height,
                                    const double diameter,
                                    const INSHierarchyIntegrator* fluid_solver,
                                    Pointer<PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Destructor.
     */
    virtual ~hagen_poiseuille_FeedbackForcer();

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
                        Pointer<PatchLevel<NDIM> > patch_level = nullptr);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    hagen_poiseuille_FeedbackForcer();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    hagen_poiseuille_FeedbackForcer(const hagen_poiseuille_FeedbackForcer& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    hagen_poiseuille_FeedbackForcer& operator=(const hagen_poiseuille_FeedbackForcer& that);

    const double d_H;
    const double d_D;
    const INSHierarchyIntegrator* const d_fluid_solver;
    Pointer<PatchHierarchy<NDIM> > d_patch_hierarchy;
};

#endif // #ifndef included_hagen_poiseuille_FeedbackForcer
