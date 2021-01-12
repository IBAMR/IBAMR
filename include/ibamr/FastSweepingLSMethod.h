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

#ifndef included_IBAMR_FastSweepingLSMethod
#define included_IBAMR_FastSweepingLSMethod

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
 * \brief Class FastSweepingLSMethod provides a fast-sweeping algorithm implementation
 * of the level set method. Specifically, this class produces a solution to the Eikonal
 * equation \f$ |\nabla Q | = 1 \f$, which produces the signed distance away from an
 * interface.
 *
 * \note The class can also compute distance function from physical domain boundary if
 * specified through input file. In presence of a physical domain wall, the distance function
 * at a grid point is D = min(distance from interface, distance from wall location).
 *
 * References
 * Zhao, H., <A HREF="http://www.ams.org/journals/mcom/2005-74-250/S0025-5718-04-01678-3/">
 * A Fast Sweeping Method For Eikonal Equations</A>
 */
class FastSweepingLSMethod : public IBAMR::LSInitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    FastSweepingLSMethod(std::string object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db = NULL,
                         bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~FastSweepingLSMethod() = default;

    /*!
     * \brief Initialize level set data using the fast-sweeping method.
     */
    void initializeLSData(int D_idx,
                          SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hierarchy_math_ops,
                          int integrator_step,
                          double time,
                          bool initial_time) override;

protected:
    // Algorithm parameters.
    bool d_consider_phys_bdry_wall = false;
    int d_wall_location_idx[2 * NDIM];

private:
    /*!
     * \brief Do one fast sweep over the hierarchy.
     */
    void fastSweep(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops, int dist_idx) const;

    /*!
     * \brief Do one fast sweep over a patch.
     */
    void fastSweep(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dist_data,
                   const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                   const SAMRAI::hier::Box<NDIM>& domain_box) const;

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
    FastSweepingLSMethod(const FastSweepingLSMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FastSweepingLSMethod& operator=(const FastSweepingLSMethod& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_FastSweepingLSMethod
