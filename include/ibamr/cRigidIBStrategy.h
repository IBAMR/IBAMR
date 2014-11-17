// Filename: cRigidIBStrategy.h
// Created on 8 Nov 2014 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

#ifndef included_cRigidIBStrategy
#define included_cRigidIBStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "ibamr/IBStrategy.h"
#include "ibamr/RigidBodyStrategy.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace xfer
{
template<int dim>
class RefineSchedule;
template<int dim>
class CoarsenSchedule;
}// namespace xfer
}// namespace SAMRAI
namespace IBTK
{
class RobinPhysBdryPatchStrategy;
}// namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
	
/*!
 * \brief Class cRigidIBStrategy is an abstract IBstrategy class which provides
 * additional support for constraint based IB methods for rigid bodies.
 */
class cRigidIBStrategy : public IBStrategy, public RigidBodyStrategy
{
//////////////////////////////////////////////////////////////////////////////
public:
	/*!
	 * \brief Constructor of the class.
	 */
	cRigidIBStrategy(
		const unsigned int parts);

	/*!
     * \brief Destructor of the class.
     */
	~cRigidIBStrategy();

	/*!
	 * \brief Spread the constraint Lagrangian force for all parts in the PetscMultiVec
	 * \f$L\f$ to the Cartesian grid at the specified time within the current
	 * time interval.
	 *
	 * \param scale Scales the Lagrangian vector before spreading.
	 */
	virtual void
	spreadForce(
		int f_data_idx,
		Vec L,
		IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
			f_prolongation_scheds,
		double data_time,
		double scale = 1.0) = 0;
	
	/*!
	 * \brief Spread the constraint Lagrangian force for a specific part in the Vec
	 * \f$L\f$ to the Cartesian grid at the specified time within the current
	 * time interval.
	 *
	 * \param scale Scales the Lagrangian vector before spreading.
	 */
	virtual void
	spreadForce(
		int f_data_idx,
		const unsigned int part,
		Vec L,
		IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
	        f_prolongation_scheds,
		double data_time,
		double scale = 1.0) = 0;
	
	/*!
	 * \brief Interpolate the Eulerian velocity to the curvilinear mesh for all parts
	 * in the PetscMultiVec \f$V\f$ at the specified time within the current 
	 * time interval.
	 *
	 * \param scale Scales the Lagrangian vector after interpolating from Eulerian grid.
	 */
	virtual void
	interpolateVelocity(
	    int u_data_idx,
		Vec V,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
			u_synch_scheds,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
			u_ghost_fill_scheds,
		double data_time,
		double scale = 1.0) = 0;
	
	/*!
	 * \brief Interpolate the Eulerian velocity to the curvilinear mesh for a specific part
	 * in the Vec \f$V\f$ at the specified time within the current time interval.
	 *
	 * \param scale Scales the Lagrangian vector after interpolating from Eulerian grid.
	 */
	virtual void
	interpolateVelocity(
		int u_data_idx,
		const unsigned int part,
		Vec V,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
			u_synch_scheds,
		const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
			u_ghost_fill_scheds,
		double data_time,
		double scale = 1.0) = 0;

//////////////////////////////////////////////////////////////////////////////
private:
	/*!
	 * \brief Copy constructor.
	 *
	 * \note This constructor is not implemented and should not be used.
	 *
	 * \param from The value to copy to this object.
	 */
	cRigidIBStrategy(
	    const cRigidIBStrategy& from);
	
	/*!
	 * \brief Assignment operator.
	 *
	 * \note This operator is not implemented and should not be used.
	 *
	 * \param that The value to assign to this object.
	 *
	 * \return A reference to this object.
	 */
	cRigidIBStrategy& operator=(
	    const cRigidIBStrategy& that);
	
};// cRigidIBStrategy
	
}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //ifndef included_cRigidIBStrategy