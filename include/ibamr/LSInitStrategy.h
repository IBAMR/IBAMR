// Filename: LSInitStrategy.h
// Created on 27 Sep 2017 by Amneet Bhalla and Nishant Nangia
//
// Copyright (c) 2002-2017, Amneet Bhalla and Nishant Nangia
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

#ifndef included_IBAMR_LSInitStrategy
#define included_IBAMR_LSInitStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "ibamr/ibamr_enums.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

namespace IBTK
{
class HierarchyMathOps;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LSInitStrategy provides a generic interface for initializing the
 * implementation details of a particular version of the level set method.
 */
class LSInitStrategy : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Constructor.
     */
    LSInitStrategy(const std::string& object_name, bool register_for_restart = true);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~LSInitStrategy();

    /*!
     * \brief Register physical boundary condition for the level set.
     */
    virtual void registerPhysicalBoundaryCondition(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* robin_bc_coef);

    /*!
     * \brief Function specifying distance function near an interface.
     */
    typedef void (*LocateInterfaceNeighborhoodFcnPtr)(int D_idx,
                                                      SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                      double time,
                                                      bool initial_time,
                                                      void* ctx);

    /*!
     * \brief Register interface neighborhood locating functions.
     */
    virtual void registerInterfaceNeighborhoodLocatingFcn(LocateInterfaceNeighborhoodFcnPtr callback, void* ctx);

    /*!
     * \brief Initialize level set data.
     */
    virtual void initializeLSData(int dst_idx,
                                  SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hierarchy_math_ops,
                                  int integrator_step,
                                  double time,
                                  bool initial_time) = 0;

    /*!
     * \brief Indicate that the class should reinitialize level set.
     */
    virtual void setReinitializeLSData(bool reinit_ls_data);

    /*!
     * Write out object state to the given database.
     *
     * An empty default implementation is provided.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    // Book-keeping.
    std::string d_object_name;
    bool d_registered_for_restart;

    // Level set order.
    LevelSetOrder d_ls_order;

    // Algorithm parameters.
    double d_abs_tol;
    int d_max_its;
    bool d_enable_logging;
    bool d_reinitialize_ls;
    int d_reinit_interval;

    // Boundary condition object for level set.
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_bc_coef;

    // Neighborhood locating functions.
    std::vector<LocateInterfaceNeighborhoodFcnPtr> d_locate_interface_fcns;
    std::vector<void*> d_locate_interface_fcns_ctx;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LSInitStrategy(const LSInitStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LSInitStrategy& operator=(const LSInitStrategy& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_LSInitStrategy
