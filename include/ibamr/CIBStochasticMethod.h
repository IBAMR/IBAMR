// Filename: CIBStochasticMethod.h
// Created on 21 Apr 2015 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith.
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
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_IBAMR_CIBStochasticMethod
#define included_IBAMR_CIBStochasticMethod

/////////////////////////////// INCLUDES /////////////////////////////////////


#include "RobinBcCoefStrategy.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/IBMethod.h"
#include "ibamr/CIBMethod.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"

namespace mu
{
class Parser;
} // namespace mu

namespace IBTK
{
class HierarchyMathsOps;
} // namespace IBTK
namespace IBAMR
{
class CIBStandardInitializer;
class CIBStaggeredStokesOperator;
} // namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CIBFEMethod is a concrete CIBStrategy and IBMethod
 * class which implements the motion of rigid bodies using the constraint
 * formulation. The immersed structure is discretized using standard IB
 * markers.
 */

class CIBStochasticMethod : public IBAMR::CIBMethod
{
public:
    /*!
     * \brief Constructor of the class.
     */
    CIBStochasticMethod(const std::string& object_name,
              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
              const int no_structures = 1,
              bool register_for_restart = true);

    /*!
     * \brief Destructor of the class.
     */
    ~CIBStochasticMethod();
    
    /*!
     * \brief setter for KT.
     */
    void setkT(double kT);
    
    /*!
     * \brief setter for Length scale.
     */
    void setLScale(double scale);
    
    /*!
     * \brief setter for RFD delta.
     */
    void setRFdelta(double scale);
    
    /*!
     * \brief getter for KT.
     */
    double getkT() const;
    
    /*!
     * \brief getter for Length scale.
     */
    double getLScale() const;
 
    /*!
     * \brief getter for RFD delta.
     */
    double getRFdelta() const;
    
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles);
    
    /*!
     * \brief Advance the positions of the Lagrangian structures
     */
    void midpointStep(double current_time, double new_time);
    
    void checkLagUpdate();
    
    /*!
     * \brief Advance the positions of the Lagrangian structure by a small increment
     */
    void moveLagrangianData(double delta);
       
    /*!
     * \brief helper function for RFD
     */
    void computeRFDforcesAndDisplacements();
    
    /*!
     * \brief sets rigid body velocity of structures at the half-time positions 
     */
    void setHalfTimeVelocity(Vec U);
    
    /*!
     * \brief resetsets rigid body velocity of structures at the half-time positions to those at current time
     */
    void resetRFDVelocity();
    
    /*!
     * \brief configuration getter
     */
    void getConfiguration(const unsigned int part, Eigen::Vector3d& C, Eigen::Quaterniond& Q);
    
     /*!
     * \brief initial configuration getter
     */
    void getInitialCOM(const unsigned int part, Eigen::Vector3d& C);
    
    bool getReject();
    
    int getNumReject();

    //////////////////////////////////////////////////////////////////////////////

protected:
    //////////////////////////////////////////////////////////////////////////////

private:
    // Boltzman Constant
    double d_kT, d_L_scale, d_rf_delta;
    bool d_reject;
    int d_num_reject;
    

}; // CIBStochasticMethod
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_CIBStochasticMethod
