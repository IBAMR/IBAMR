// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_BrinkmanAdvDiffSemiImplicitHierarchyIntegrator
#define included_IBAMR_BrinkmanAdvDiffSemiImplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/config.h"

#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"

#include "HierarchyFaceDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <set>
#include <string>

namespace IBAMR
{
class BrinkmanAdvDiffBcHelper;
class ConvectiveOperator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class BrinkmanAdvDiffSemiImplicitHierarchyIntegrator manages the spatial
 * discretization and time integration of scalar- and vector-valued quantities
 * whose dynamics are governed by the advection-diffusion equation. If the variable
 * is registered with penalization boundary conditions, this class computes the additional
 * diffusion coefficient and the forcing term from BrinkmanAdvDiffBcHelper class
 * and solves a penalized advection-diffusion equation for that variable.
 *
 * Each quantity \f$ Q \f$ managed by the integrator may have a unique diffusion
 * coefficient \f$ \kappa \f$ and damping coefficient \f$ \lambda \f$, and may
 * optionally have a forcing term \f$ F \f$.  Additionally, a different
 * advection velocity may be used with each quantity registered with the
 * integrator.
 *
 * This hierarchy integrator advances all levels of the patch hierarchy
 * synchronously in time.  In particular, subcycling in time is \em not
 * performed.
 *
 * Various options are available for the spatial and temporal discretizations.
 *
 * \see HierarchyIntegrator
 * \see SAMRAI::mesh::StandardTagAndInitStrategy
 * \see SAMRAI::algs::TimeRefinementIntegrator
 * \see SAMRAI::algs::TimeRefinementLevelStrategy
 */
class BrinkmanAdvDiffSemiImplicitHierarchyIntegrator : public AdvDiffSemiImplicitHierarchyIntegrator
{
public:
    /*!
     * The constructor for class BrinkmanAdvDiffSemiImplicitHierarchyIntegrator sets
     * some default values, reads in configuration information from input and
     * restart databases, and registers the integrator object with the restart
     * manager when requested.
     */
    BrinkmanAdvDiffSemiImplicitHierarchyIntegrator(const std::string& object_name,
                                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                   bool register_for_restart = true);

    /*!
     * Empty destructor.
     */
    ~BrinkmanAdvDiffSemiImplicitHierarchyIntegrator() = default;

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;
    /*!
     * \brief Register BrinkmanAdvDiffBcHelper object to add Brinkman penalization terms to the advection-diffusion
     * solver.
     */
    void registerBrinkmanAdvDiffBcHelper(SAMRAI::tbox::Pointer<IBAMR::BrinkmanAdvDiffBcHelper> brinkman_penalization);

    /*!
     * \brief Get the BrinkmanAdvDiffBcHelper object registered with this class.
     */
    const SAMRAI::tbox::Pointer<IBAMR::BrinkmanAdvDiffBcHelper>& getBrinkmanPenalization() const
    {
        return d_brinkman_penalization;
    } // getBrinkmanPenalization

protected:
    /*!
     * Additional variables required for Brinkman penalization
     */
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > >
        d_Q_Cb_map, d_Q_Cb_rhs_map, d_Q_Fb_map;
    std::map<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> >,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > >
        d_Q_Db_map, d_Q_Db_rhs_map;

    /*!
     * Flag to zero out the temporal term contribution when the Brinkman approach
     * is used.
     */
    bool d_brinkman_time_independent = false;

    /*!
     * Brinkman penalization object registred with this integrator.
     */
    SAMRAI::tbox::Pointer<IBAMR::BrinkmanAdvDiffBcHelper> d_brinkman_penalization;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    BrinkmanAdvDiffSemiImplicitHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BrinkmanAdvDiffSemiImplicitHierarchyIntegrator(const BrinkmanAdvDiffSemiImplicitHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BrinkmanAdvDiffSemiImplicitHierarchyIntegrator&
    operator=(const BrinkmanAdvDiffSemiImplicitHierarchyIntegrator& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_BrinkmanAdvDiffSemiImplicitHierarchyIntegrator
