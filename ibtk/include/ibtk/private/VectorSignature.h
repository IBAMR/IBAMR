// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#ifndef included_IBTK_VectorSignature
#define included_IBTK_VectorSignature

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <tbox/Pointer.h>

#include <IntVector.h>
#include <PatchHierarchy.h>
#include <SAMRAIVectorReal.h>
#include <Variable.h>

#include <cstddef>
#include <typeindex>
#include <typeinfo>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
namespace detail
{
/*!
 * \brief The structure of a SAMRAIVectorReal: its hierarchy, level range, and, for each component, its variable, patch
 * data factory type, depth, ghost width, control-volume index, and descriptor index. Solvers that reuse state built
 * from a vector compare the signatures of successive vectors. The hierarchy and the variables of a stored signature
 * stay alive until the signature is reset.
 */
class VectorSignature
{
public:
    /*!
     * \brief Construct an empty signature, which differs from the signature of every vector.
     */
    VectorSignature() = default;

    /*!
     * \brief Construct the signature of \a vec.
     */
    explicit VectorSignature(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& vec);

    /*!
     * \brief Return true if \a other has the same hierarchy, level range, and components, including the descriptor
     * indices of the data.
     */
    bool operator==(const VectorSignature& other) const;

    /*!
     * \brief Return true if \a other has the same hierarchy and components, ignoring the descriptor indices of the
     * data and the level range. A clone of a vector with this signature has the structure of a clone of \a other.
     */
    bool hasSameStructure(const VectorSignature& other) const;

private:
    /*!
     * \brief The signature of one component.
     */
    struct Component
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> variable;
        std::type_index factory_type = typeid(void);
        std::size_t memory_size =
            0; // the memory of a one-cell box, which gives the depth of the data given the type and ghost width
        SAMRAI::hier::IntVector<NDIM> ghost_width;
        int control_volume_idx = -1;
        int data_idx = -1;
    };

    bool hasSameStructure(const VectorSignature& other, bool compare_data_indices) const;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    int d_coarsest_ln = -1;
    int d_finest_ln = -1;
    std::vector<Component> d_components;
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_VectorSignature
