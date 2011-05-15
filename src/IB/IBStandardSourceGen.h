// Filename: IBStandardSourceGen.h
// Created on 28 Apr 2011 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
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

#ifndef included_IBStandardSourceGen
#define included_IBStandardSourceGen

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/IBLagrangianSourceStrategy.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBStandardSourceGen provides support for distributed internal
 * fluid sources/sinks.
 */
class IBStandardSourceGen
    : public virtual IBLagrangianSourceStrategy
{
public:
    /*!
     * \brief Default constructor.
     */
    IBStandardSourceGen();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBStandardSourceGen();

    /*!
     * \brief Return a reference to the vector of source strengths.
     *
     * \note Users \em must \em not change the size of this vector.
     */
    std::vector<double>&
    getSourceStrengths(
        const int ln);

    /*!
     * \brief Return a const reference to the vector of source strengths.
     */
    const std::vector<double>&
    getSourceStrengths(
        const int ln) const;

    /*!
     * \brief Return a const reference to the vector of source pressures.
     */
    const std::vector<double>&
    getSourcePressures(
        const int ln) const;

    /*!
     * \brief Return the number of sources/sinks.
     */
    int
    getNumSources(
        const int ln) const;

    /*!
     * \brief Return a reference to the vector of source names.
     *
     * \note Users \em must \em not change the size of this vector.
     */
    std::vector<std::string>&
    getSourceNames(
        const int ln);

    /*!
     * \brief Return a const reference to the vector of source names.
     */
    const std::vector<std::string>&
    getSourceNames(
        const int ln) const;

    /*!
     * \brief Return a reference to the vector of source radii.
     *
     * \note Users \em must \em not change the size of this vector.
     */
    std::vector<double>&
    getSourceRadii(
        const int ln);

    /*!
     * \brief Return a const reference to the vector of source radii.
     */
    const std::vector<double>&
    getSourceRadii(
        const int ln) const;

    /*!
     * \brief Setup the data needed to compute source/sink data on the specified
     * level of the patch hierarchy.
     *
     * \note A default empty implementation is provided.
     */
    virtual void
    initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool initial_time,
        IBTK::LDataManager* const l_data_manager);

    /*!
     * \brief Specify the number of distributed internal sources or sinks.
     *
     * \note The return value must be the \em total number of internal
     * sources/sinks in the \em entire computational domain.  This implies that
     * the return value must be \em identical on each MPI process.
     */
    virtual int
    getNumSources(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        IBTK::LDataManager* const l_data_manager);

    /*!
     * \brief Compute the source locations for each of the distributed internal
     * sources or sinks.
     *
     * \note Implementations of this method \em must compute the same values for
     * \a X_src on \em each MPI process.  That is to say, \a X_src must provide
     * the location of all of the distributed sources/sinks.
     */
    virtual void
    getSourceLocations(
        std::vector<std::vector<double> >& X_src,
        std::vector<double>& r_src,
        SAMRAI::tbox::Pointer<IBTK::LData> X_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        IBTK::LDataManager* const l_data_manager);

    /*!
     * \brief Set the normalized pressures at the sources.
     */
    virtual void
    setSourcePressures(
        const std::vector<double>& P_src,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        IBTK::LDataManager* const l_data_manager);

    /*!
     * \brief Compute the source strengths for each of the distributed internal
     * sources or sinks.
     *
     * \note Implementations of this method \em must compute the same values for
     * \a Q_src on \em each MPI process.  That is to say, \a Q_src must provide
     * the strengths of all of the distributed sources/sinks.
     */
    virtual void
    computeSourceStrengths(
        std::vector<double>& Q_src,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double data_time,
        IBTK::LDataManager* const l_data_manager);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBStandardSourceGen(
        const IBStandardSourceGen& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBStandardSourceGen&
    operator=(
        const IBStandardSourceGen& that);

    /*
     * Source/sink data.
     */
    std::vector<int> d_n_src;
    std::vector<std::vector<std::string> >d_source_names;
    std::vector<std::vector<double> > d_r_src;
    std::vector<std::vector<int> > d_num_perimeter_nodes;
    std::vector<std::vector<double> > d_Q_src, d_P_src;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBStandardSourceGen.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBStandardSourceGen
