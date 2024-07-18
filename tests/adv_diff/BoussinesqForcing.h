// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_BoussinesqForcing
#define included_BoussinesqForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////
/*!
 * \brief Class BoussinesqForcing provides forcing for the momentum equations
 * based on the Boussinesq approximation.
 */
class BoussinesqForcing : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    BoussinesqForcing(SAMRAIPointer<SAMRAI::hier::VariableNd> T_var,
                      SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
                      IBTK::SAMRAIPointer<CellVariableNd<double> > ls_inner_solid_var,
                      IBTK::SAMRAIPointer<CellVariableNd<double> > ls_outer_solid_var,
                      SAMRAIPointer<Database> input_db);

    /*!
     * \brief Empty destructor.
     */
    ~BoussinesqForcing();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete BoussinesqForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 SAMRAIPointer<SAMRAI::hier::VariableNd> var,
                                 SAMRAIPointer<PatchHierarchyNd> hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAIPointer<SAMRAI::hier::VariableNd> var,
                        SAMRAIPointer<PatchNd> patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAIPointer<PatchLevelNd> patch_level = SAMRAIPointer<PatchLevelNd>(nullptr));

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BoussinesqForcing(const BoussinesqForcing& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BoussinesqForcing& operator=(const BoussinesqForcing& that) = delete;

    /*
     * Pointer to the non-dimensional temperature.
     */
    SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > d_T_var;

    /*
     * Pointer to the adv-diff solver.
     */
    SAMRAIPointer<AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;

    /*
     * Vector storing the pointer to level set variables.
     */
    std::vector<IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > > d_ls_solid_vars;

    /*
     * A variable and index to store the sum of mask functions i.e., \f$ \sum \chi_i\f$
     */
    IBTK::SAMRAIPointer<SAMRAI::pdat::SideVariableNd<double> > d_chi_var;
    int d_chi_idx = IBTK::invalid_index;

    /*
     * Non-dimensional numbers
     */
    double d_rayleigh_number, d_prandtl_number;
};

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_BoussinesqForcing
