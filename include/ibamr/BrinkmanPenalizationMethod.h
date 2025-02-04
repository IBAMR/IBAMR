// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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

#ifndef included_BrinkmanPenalizationMethod
#define included_BrinkmanPenalizationMethod

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/BrinkmanPenalizationStrategy.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
class HierarchyIntegrator;
} // namespace IBTK
namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief BrinkmanPenalizationMethod is a simplified (and a lightweight) implementation of the
 * BrinkmanPenalizationStrategy class that does not require any specific HierarchyIntegrator.
 */
class BrinkmanPenalizationMethod : public BrinkmanPenalizationStrategy
{
public:
    /*
     * \brief Constructor of the class.
     */
    BrinkmanPenalizationMethod(std::string object_name,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> time_integrator,
                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                               bool register_for_restart = true);

    /*
     * \brief Destructor of the class.
     */
    ~BrinkmanPenalizationMethod() = default;

    /*
     * \brief Register the level set function of the solid body.
     */
    void registerSolidLevelSet(int ls_current_idx,
                               int ls_new_idx,
                               int ls_scratch_idx,
                               SAMRAI::solv::RobinBcCoefStrategy<NDIM>* ls_bc_coef);

    /*
     * \brief Get the current patch index of the solid level set function.
     */
    int getLevelSetCurrentPatchIndex()
    {
        return d_ls_current_idx;
    } // getLevelSetCurrentPatchIndex

    /*
     * \brief Get the new patch index of the solid level set function.
     */
    int getLevelSetNewPatchIndex()
    {
        return d_ls_new_idx;
    } // getLevelSetNewPatchIndex

    /*
     * \brief Get the scratch patch index of the solid level set function.
     */
    int getLevelSetScratchPatchIndex()
    {
        return d_ls_scratch_idx;
    } // getLevelSetScratchPatchIndex

    /*!
     * \brief Compute the desired rigid body velocity in the Brinkman penalized (solid) zone.
     */
    void computeBrinkmanVelocity(int u_idx, double time, int cycle_num) override;

    /*!
     * \brief Demarcate the Brinkman zone.
     */
    void demarcateBrinkmanZone(int u_idx, double time, int cycle_num) override;

    /*!
     * \brief Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*!
     * \brief Pointer to the time integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> d_time_integrator;

    // Number of interface cells to compute the Heaviside function
    double d_num_interface_cells = 2.0;

    /*!
     * \brief Patch index and boundary condition object for the solid level set variable.
     */
    int d_ls_current_idx = IBTK::invalid_index, d_ls_new_idx = IBTK::invalid_index,
        d_ls_scratch_idx = IBTK::invalid_index;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_ls_bc_coef = nullptr;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BrinkmanPenalizationMethod(const BrinkmanPenalizationMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BrinkmanPenalizationMethod& operator=(const BrinkmanPenalizationMethod& that) = delete;

    /*!
     * \brief Get options from input database.
     */
    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_BrinkmanPenalizationMethod
