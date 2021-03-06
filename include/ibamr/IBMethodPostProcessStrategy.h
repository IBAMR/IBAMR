// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_IBMethodPostProcessStrategy
#define included_IBAMR_IBMethodPostProcessStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "tbox/DescribedClass.h"

#include <vector>

namespace IBTK
{
class LData;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
} // namespace SAMRAI
/////////////////////////////// INCLUDES /////////////////////////////////////

namespace IBAMR
{
class IBMethod;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBMethodPostProcessStrategy provides a generic interface for
 * specifying post-processing code for use in an IB computation.
 */
class IBMethodPostProcessStrategy : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBMethodPostProcessStrategy() = default;

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBMethodPostProcessStrategy() = default;

    /*!
     * \brief Post-process data on the patch hierarchy.
     */
    virtual void postprocessData(int u_idx,
                                 int p_idx,
                                 int f_idx,
                                 const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& F_data,
                                 const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& X_data,
                                 const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& U_data,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 int coarsest_level_number,
                                 int finest_level_number,
                                 double data_time,
                                 IBMethod* ib_method) = 0;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMethodPostProcessStrategy(const IBMethodPostProcessStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMethodPostProcessStrategy& operator=(const IBMethodPostProcessStrategy& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBMethodPostProcessStrategy
