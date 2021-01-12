// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_RelaxationLSMethod
#define included_IBAMR_RelaxationLSMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/LSInitStrategy.h"
#include "ibamr/ibamr_enums.h"

#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace pdat
{
template <int DIM, class TYPE>
class CellData;
} // namespace pdat
namespace hier
{
template <int DIM>
class Box;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class RelaxationLSMethod provides a relaxation algorithm implementation
 * of the level set method. Specifically, this class iterates (to steady-state) the PDE
 * \f$\frac{\partial Q}{\partial \tau}+sgn(Q^0)(|\nabla Q | - 1) = 0\f$, which produces a solution
 * to the Eikonal equation \f$ |\nabla Q | = 1 \f$. The solution of the Eikonal equation
 * produces the signed distance away from an interface.
 *
 * An optional mass constraint can be applied to certain high-order LS schemes by
 * first computing an intermediate \f$\tilde{Q}\f$ and projecting in the following way:

 * \f$Q = \tilde{Q} + \lambda_{ij}H'(Q^0)|\nabla Q^0|\f$.

 * This greatly improves the volume conservation of the interface while still
 * retaining the order of accuracy and signed distance property. The mass conservation
 * constraint assumes that \f$Q^0\f$ is already close to a signed distance function and
 * is hence, by default, disabled at initial time.
 *
 *
 * References
 * Min, C., <A HREF="http://www.sciencedirect.com/science/article/pii/S0021999109007189">
 * On reinitializing level set functions</A>
 *
 * Sussman, M. and Fatemi, E., <A HREF="http://epubs.siam.org/doi/abs/10.1137/S1064827596298245">
 * An Efficient, Interface-Preserving Level Set Redistancing Algorithm and Its Application to Interfacial Incompressible
 * Fluid Flow </A>
 */
class RelaxationLSMethod : public IBAMR::LSInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    RelaxationLSMethod(std::string object_name,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db = NULL,
                       bool register_for_restart = true);

    /*!
     * \name Implementation of IBAMR::LSInitStrategy interface.
     */
    //\{

    /*!
     * \brief Initialize level set data using the relaxation method.
     */
    void initializeLSData(int D_idx,
                          SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hierarchy_math_ops,
                          int integrator_step,
                          double time,
                          bool initial_time) override;

    //\}

    /*!
     * \brief Indicate that the class should apply the mass constraint.
     */
    void setApplyMassConstraint(bool apply_mass_constraint);

    /*!
     * \brief Indicate that the class should apply the subcell fix.
     */
    void setApplySubcellFix(bool apply_subcell_fix);

    /*!
     * \brief Indicate that the class should apply the sign fix near interface points.
     */
    void setApplySignFix(bool apply_sign_fix);

    /*!
     * \brief Hard set the ghost cell width for the level set variable.
     */
    void setLSGhostCellWidth(int D_gcw);

    /*!
     * \brief Indicate that the class should apply the volume shift.
     */
    void setApplyVolumeShift(bool apply_volume_shift);

protected:
    // Flag for applying the mass constraint
    bool d_apply_mass_constraint = false;

    // Flag for applying subcell fix
    bool d_apply_subcell_fix = false;

    // Flag for applying sign fix near interface points
    bool d_apply_sign_fix = false;

    // Ghost cell width for level set variable
    int d_D_gcw = -1;

    // Flag for applying the volume shift
    bool d_apply_volume_shift = false;

    // Initial volume of the level set domain
    double d_init_ls_vol;

    // Relaxation weight parameter
    double d_alpha = 1.0;

private:
    /*!
     * \brief Do one relaxation step over the hierarchy.
     */
    void relax(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
               int dist_idx,
               int dist_init_idx,
               const int iter) const;

    /*!
     * \brief Do one relaxation step over a patch.
     */
    void relax(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_data,
               const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_init_data,
               const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
               const int iter) const;

    /*!
     * \brief Compute the Hamiltonian of the indicator field over the hierarchy
     */
    void computeInitialHamiltonian(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                   int ham_init_idx,
                                   int dist_init_idx) const;

    /*!
     * \brief Compute the hamiltonian of the indicator field field over a patch
     */
    void computeInitialHamiltonian(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > ham_init_data,
                                   const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_init_data,
                                   const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Apply the mass constraint over the hierarchy
     */
    void applyMassConstraint(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             int dist_idx,
                             int dist_copy_idx,
                             int dist_init_idx,
                             int ham_init_idx) const;

    /*!
     * \brief Apply the mass constraint over a patch
     */
    void applyMassConstraint(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_data,
                             const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_copy_data,
                             const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_init_data,
                             const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > ham_init_data,
                             const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Compute the volume of a region demarcated by a level set variable
     */
    double
    computeRegionVolume(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops, int hs_phi_idx, int phi_idx) const;

    /*!
     * \brief Apply the volume shift over the hierarchy
     */
    void applyVolumeShift(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                          int dist_idx,
                          int dist_copy_idx,
                          double dV) const;

    /*!
     * \brief Apply the volume shift over a patch
     */
    void applyVolumeShift(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_data,
                          const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_copy_data,
                          const double dV,
                          const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    RelaxationLSMethod(const RelaxationLSMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    RelaxationLSMethod& operator=(const RelaxationLSMethod& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_RelaxationLSMethod
