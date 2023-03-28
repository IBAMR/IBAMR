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

#ifndef included_IBTK_LMesh
#define included_IBTK_LMesh

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/DescribedClass.h"

#include <string>
#include <vector>

namespace IBTK
{
class LNode;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LMesh is a collection of LNode objects.
 */
class LMesh : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    LMesh(std::string object_name, std::vector<LNode*> local_nodes, std::vector<LNode*> ghost_nodes);

    /*!
     * \brief Destructor.
     */
    virtual ~LMesh() = default;

    /*!
     * \brief Return a const reference to the set of local LNode objects.
     */
    const std::vector<LNode*>& getLocalNodes() const;

    /*!
     * \brief Return a const reference to the set of local ghost LNode objects.
     */
    const std::vector<LNode*>& getGhostNodes() const;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LMesh(const LMesh& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LMesh& operator=(const LMesh& that) = delete;

    const std::string& d_object_name;
    const std::vector<LNode*> d_local_nodes;
    const std::vector<LNode*> d_ghost_nodes;
};

} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LMesh-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LMesh
