// Filename: IBLagrangianSourceStrategy.h
// Created on 18 Jun 2005 by Boyce Griffith
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

#ifndef included_IBLagrangianSourceStrategy
#define included_IBLagrangianSourceStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "ibtk/ibtk_utilities.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class LData;
class LDataManager;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBLagrangianForceStrategy provides a generic interface for
 * specifying the positions and magnitudes of distributed internal fluid
 * source-sinks.
 */
class IBLagrangianSourceStrategy : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    IBLagrangianSourceStrategy();

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBLagrangianSourceStrategy();

    /*!
     * \brief Set the current and new times for the present timestep.
     *
     * \note A default empty implementation is provided.
     */
    virtual void setTimeInterval(double current_time, double new_time);

    /*!
     * \brief Setup the data needed to compute source/sink data on the specified
     * level of the patch hierarchy.
     *
     * \note A default empty implementation is provided.
     */
    virtual void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     int level_number,
                                     double init_data_time,
                                     bool initial_time,
                                     IBTK::LDataManager* l_data_manager);

    /*!
     * \brief Specify the number of distributed internal sources or sinks.
     *
     * \note The return value must be the \em total number of internal
     * sources/sinks in the \em entire computational domain.  This implies that
     * the return value must be \em identical on each MPI process.
     */
    virtual unsigned int getNumSources(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       int level_number,
                                       double data_time,
                                       IBTK::LDataManager* l_data_manager) = 0;

    /*!
     * \brief Compute the source locations for each of the distributed internal
     * sources or sinks.
     *
     * \note Implementations of this method \em must compute the same values for
     * \a X_src on \em each MPI process.  That is to say, \a X_src must provide
     * the location of all of the distributed sources/sinks.
     */
    virtual void getSourceLocations(std::vector<IBTK::Point>& X_src,
                                    std::vector<double>& r_src,
                                    SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                    int level_number,
                                    double data_time,
                                    IBTK::LDataManager* l_data_manager) = 0;

    /*!
     * \brief Set the normalized pressures at the sources.
     */
    virtual void setSourcePressures(const std::vector<double>& P_src,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                    int level_number,
                                    double data_time,
                                    IBTK::LDataManager* l_data_manager) = 0;

    /*!
     * \brief Compute the source strengths for each of the distributed internal
     * sources or sinks.
     *
     * \note Implementations of this method \em must compute the same values for
     * \a Q_src on \em each MPI process.  That is to say, \a Q_src must provide
     * the strengths of all of the distributed sources/sinks.
     */
    virtual void computeSourceStrengths(std::vector<double>& Q_src,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                        int level_number,
                                        double data_time,
                                        IBTK::LDataManager* l_data_manager) = 0;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBLagrangianSourceStrategy(const IBLagrangianSourceStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBLagrangianSourceStrategy& operator=(const IBLagrangianSourceStrategy& that);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBLagrangianSourceStrategy
