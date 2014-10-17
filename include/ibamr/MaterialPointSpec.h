// Filename: MaterialPointSpec.h
// Created on 16 Oct 2012 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#ifndef included_MaterialPointSpec
#define included_MaterialPointSpec

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <vector>

#include "ibtk/Streamable.h"
#include "ibtk/StreamableFactory.h"
#include "libmesh/id_types.h"
#include "tbox/Pointer.h"

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
    int getStreamableClassID() const;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    size_t getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    void packStream(SAMRAI::tbox::AbstractStream& stream);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    MaterialPointSpec(const MaterialPointSpec& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    MaterialPointSpec& operator=(const MaterialPointSpec& that);

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
        ~Factory();

        /*!
         * \brief Return the unique identifier used to specify the
         * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
         * extract MaterialPointSpec objects from data streams.
         */
        int getStreamableClassID() const;

        /*!
         * \brief Set the unique identifier used to specify the
         * IBTK::StreamableFactory object used by the IBTK::StreamableManager to
         * extract MaterialPointSpec objects from data streams.
         */
        void setStreamableClassID(int class_id);

        /*!
         * \brief Build an MaterialPointSpec object by unpacking data from the
         * data stream.
         */
        SAMRAI::tbox::Pointer<IBTK::Streamable> unpackStream(SAMRAI::tbox::AbstractStream& stream,
                                                             const SAMRAI::hier::IntVector<NDIM>& offset);

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
        Factory(const Factory& from);

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        Factory& operator=(const Factory& that);

        friend class MaterialPointSpec;
    };
    typedef MaterialPointSpec::Factory MaterialPointSpecFactory;
};
} // namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibamr/private/MaterialPointSpec-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_MaterialPointSpec
