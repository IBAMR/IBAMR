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

#ifndef included_IBAMR_LSInitStrategy
#define included_IBAMR_LSInitStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/ibamr_enums.h"

#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include <string>
#include <vector>

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
    LSInitStrategy(std::string object_name, bool register_for_restart = true);

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
    using LocateInterfaceNeighborhoodFcnPtr = void (*)(int D_idx,
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
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    // Book-keeping.
    std::string d_object_name;
    bool d_registered_for_restart;

    // Level set order.
    LevelSetOrder d_ls_order = FIRST_ORDER_LS;

    // Algorithm parameters.
    double d_abs_tol = 1.0e-5;
    int d_max_its = 100;
    bool d_enable_logging = false;
    bool d_reinitialize_ls = false;
    int d_reinit_interval = 0;

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
    LSInitStrategy(const LSInitStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LSInitStrategy& operator=(const LSInitStrategy& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_LSInitStrategy
