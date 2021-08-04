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

#ifndef included_IBAMR_MaterialPointSpec
#define included_IBAMR_MaterialPointSpec

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include "ibtk/Streamable.h"
#include "ibtk/StreamableFactory.h"

#include "tbox/Pointer.h"

#include "libmesh/id_types.h"

#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class IntVector;
} // namespace hier
namespace tbox
{
class AbstractStream;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class MaterialPointSpec encapsulates data necessary to define the
 * properties associated with an immersed material point.
 */
class MaterialPointSpec : public IBTK::Streamable
{
public:
    /*!
     * \brief Register this class and its factory class with the singleton
     * IBTK::StreamableManager object.  This method must be called before any
     * MaterialPointSpec objects are created.
     *
     * \note This method is collective on all MPI processes.  This is done to
     * ensure that all processes employ the same class ID for the
     * MaterialPointSpec class.
     */
    static void registerWithStreamableManager();

    /*!
     * \brief Returns a boolean indicating whether the class has been registered
     * with the singleton IBTK::StreamableManager object.
     */
    static bool getIsRegisteredWithStreamableManager();

    /*!
     * The unique class ID for this object type assigned by the
     * IBTK::StreamableManager.
     */
    static int STREAMABLE_CLASS_ID;

    /*!
     * \brief Default constructor.
     */
    MaterialPointSpec(int point_idx = -1,
                      double weight = 0.0,
                      libMesh::subdomain_id_type subdomain_id = 0,
                      const std::vector<double>& internal_vars = std::vector<double>());

    /*!
     * \brief Destructor.
     */
    ~MaterialPointSpec();

    /*!
     * \return A const reference to the point index.
     */
    const int& getPointIndex() const;

    /*!
     * \return A non-const reference to the point index.
     */
    int& getPointIndex();

    /*!
     * \return A const reference to the weighting factor associated with the
     * material point.
     */
    const double& getWeight() const;

    /*!
     * \return A non-const reference to the weighthting factor associated with
     * the material point.
     */
    double& getWeight();

    /*!
     * \return A const reference to the subdomain_id associated with the
     * material point.
     */
    const libMesh::subdomain_id_type& getSubdomainId() const;

    /*!
     * \return A non-const reference to the subdomain_id associated with the
     * material point.
     */
    libMesh::subdomain_id_type& getSubdomainId();

    /*!
     * \return A const reference to the internal state variables associated with
     * the material point.
     */
    const std::vector<double>& getInternalVariables() const;

    /*!
     * \return A non-const reference to the internal state variables associated
     * with the material point.
     */
    std::vector<double>& getInternalVariables();

    /*!
     * \brief Return the unique identifier used to specify the
     * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
     * extract Streamable objects from data streams.
     */
    int getStreamableClassID() const override;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t getDataStreamSize() const override;

    /*!
     * \brief Pack data into the output stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream) override;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    MaterialPointSpec(const MaterialPointSpec& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    MaterialPointSpec& operator=(const MaterialPointSpec& that) = delete;

    /*!
     * Material point data.
     */
    int d_point_idx;
    double d_weight;
    libMesh::subdomain_id_type d_subdomain_id;
    std::vector<double> d_internal_vars;

    /*!
     * \brief A factory class to rebuild MaterialPointSpec objects from
     * SAMRAI::tbox::AbstractStream data streams.
     */
    class Factory : public IBTK::StreamableFactory
    {
    public:
        /*!
         * \brief Destructor.
         */
        ~Factory() = default;

        /*!
         * \brief Return the unique identifier used to specify the
         * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
         * extract MaterialPointSpec objects from data streams.
         */
        int getStreamableClassID() const override;

        /*!
         * \brief Set the unique identifier used to specify the
         * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
         * extract MaterialPointSpec objects from data streams.
         */
        void setStreamableClassID(int class_id) override;

        /*!
         * \brief Build an MaterialPointSpec object by unpacking data from the
         * data stream.
         */
        SAMRAI::tbox::Pointer<IBTK::Streamable> unpackStream(SAMRAI::tbox::AbstractStream& stream,
                                                             const SAMRAI::hier::IntVector<NDIM>& offset) override;

    private:
        /*!
         * \brief Default constructor.
         */
        Factory();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        Factory(const Factory& from) = delete;

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        Factory& operator=(const Factory& that) = delete;

        friend class MaterialPointSpec;
    };
    using MaterialPointSpecFactory = MaterialPointSpec::Factory;
};
} // namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibamr/private/MaterialPointSpec-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBAMR_HAVE_LIBMESH
#endif //#ifndef included_IBAMR_MaterialPointSpec
