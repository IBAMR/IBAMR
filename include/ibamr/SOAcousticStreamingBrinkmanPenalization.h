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

#ifndef included_SOAcousticStreamingBrinkmanPenalization
#define included_SOAcousticStreamingBrinkmanPenalization

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/BrinkmanPenalizationMethod.h"

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
 * \brief SOAcousticStreamingBrinkmanPenalization class imposes Stokes' drift velocity in the immersed
 * domain using the Brinkman penalization technique.
 */
class SOAcousticStreamingBrinkmanPenalization : public BrinkmanPenalizationMethod
{
public:
    /*
     * \brief Constructor of the class.
     */
    SOAcousticStreamingBrinkmanPenalization(std::string object_name,
                                            SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> time_integrator,
                                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                            bool register_for_restart = true);

    /*
     * \brief Destructor of the class.
     */
    ~SOAcousticStreamingBrinkmanPenalization() = default;

    /*
     * \brief Set first-order velocity patch data indices.
     */
    void setFOVelocityPatchDataIndices(int U1_real_idx, int U1_imag_idx);

    /*
     * \brief Set acoustic angular frequency.
     */
    void setAcousticAngularFrequency(double omega);

    /*!
     * \brief Impose Stokes drift velocity in the Brinkman penalized (solid) zone.
     */
    void computeBrinkmanVelocity(int u_idx, double time, int cycle_num) override;

protected:
    /*!
     * \brief Patch indices for first order velocity.
     */
    int d_U1_real_idx = IBTK::invalid_index, d_U1_imag_idx = IBTK::invalid_index;

    /*
     * \brief Acoustic angular frequency.
     */
    double d_acoustic_freq = std::numeric_limits<double>::signaling_NaN();

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SOAcousticStreamingBrinkmanPenalization(const SOAcousticStreamingBrinkmanPenalization& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SOAcousticStreamingBrinkmanPenalization& operator=(const SOAcousticStreamingBrinkmanPenalization& that) = delete;

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

#endif // #ifndef included_SOAcousticStreamingBrinkmanPenalization
