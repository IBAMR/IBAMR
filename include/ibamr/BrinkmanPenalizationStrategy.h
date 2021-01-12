// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#ifndef included_BrinkmanPenalizationStrategy
#define included_BrinkmanPenalizationStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include <limits>
#include <string>
#include <utility>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief BrinkmanPenalizationStrategy is an abstract class that provides an interface
 * to implement Brinkman penalization body force in the momentum equation.
 */
class BrinkmanPenalizationStrategy : public SAMRAI::tbox::Serializable
{
public:
    /*
     * \brief Constructor of the class.
     */
    BrinkmanPenalizationStrategy(std::string object_name, bool register_for_restart = true);

    /*
     * \brief Destructor of the class.
     */
    virtual ~BrinkmanPenalizationStrategy();

    /*!
     * \brief Set the time interval in which Brinkman velocity is computed.
     */
    virtual void setTimeInterval(double current_time, double new_time);

    /*!
     * \brief Preprocess routine before computing Brinkman penalization related terms.
     *
     * \note An empty implementation is provided by default.
     */
    virtual void preprocessComputeBrinkmanPenalization(double current_time, double new_time, int num_cycles);

    /*!
     * \brief Compute the Brinkman penalized velocity term.
     */
    virtual void computeBrinkmanVelocity(int u_idx, double time, int cycle_num) = 0;

    /*!
     * \brief Demarcate the Brinkman penalization zone with Brinkman penalty term.
     */
    virtual void demarcateBrinkmanZone(int u_idx, double time, int cycle_num) = 0;

    /*!
     * \brief Postprocess routine after computing Brinkman penalization related terms.
     *
     * \note An empty implementation is provided by default.
     */
    virtual void postprocessComputeBrinkmanPenalization(double current_time, double new_time, int num_cycles);

    /*!
     * \brief Set Brinkman penalization penalty factor.
     */
    virtual void setBrinkmanCoefficient(double chi);

    /*!
     * \brief Write out object state to the given database.
     *
     * \note An empty default implementation is provided.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*
     * \brief Get the name of the object.
     */
    const std::string& getName() const
    {
        return d_object_name;
    } // getName

    /*
     * \brief Get the Brinkman coefficient.
     */
    double getBrinkmanCoefficient() const
    {
        return d_chi;
    } // getBrinkmanPenaltyFactor

    /*
     * \brief Get the current time interval \f$ [t^{n+1}, t^n] \f$ in which Brinkman
     * velocity is computed.
     */
    std::pair<double, double> getCurrentTimeInterval() const
    {
        return std::make_pair(d_new_time, d_current_time);
    } // getCurrentTimeInterval

protected:
    /*!
     * Book-keeping.
     */
    std::string d_object_name;

    /*
     * A boolean value indicating whether the class is registered with the
     * restart database.
     */
    bool d_registered_for_restart;

    /*
     * Time interval
     */
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN();

    /*
     * Brinkman coefficient.
     */
    double d_chi = 1e8;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BrinkmanPenalizationStrategy(const BrinkmanPenalizationStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BrinkmanPenalizationStrategy& operator=(const BrinkmanPenalizationStrategy& that) = delete;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_BrinkmanPenalizationStrategy
