// Filename: BrinkmanPenalizationStrategy.h
// Created on 04 Dec 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2018, Amneet Bhalla
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_BrinkmanPenalizationStrategy
#define included_BrinkmanPenalizationStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <limits>
#include <string>

#include "tbox/Serializable.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief BrinkmanPenalizationStrategy is an abstract class that provides an interface
 * to implement Brinkman penalization body force in the momentum equation.
 */
class BrinkmanPenalizationStrategy : public SAMRAI::tbox::Serializable
{
    ////////////////////////////// PUBLIC ////////////////////////////////////////
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

    /////////////////////////////// PROTECTED ////////////////////////////////////

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

    /////////////////////////////// PROTECTED ////////////////////////////////////

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
