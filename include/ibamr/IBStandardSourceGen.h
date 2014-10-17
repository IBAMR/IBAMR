// Filename: IBStandardSourceGen.h
// Created on 28 Apr 2011 by Boyce Griffith
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

#ifndef included_IBStandardSourceGen
#define included_IBStandardSourceGen

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "ibamr/IBLagrangianSourceStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

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
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBStandardSourceGen provides support for distributed internal
 * fluid sources/sinks.
 */
class IBStandardSourceGen : public IBLagrangianSourceStrategy, public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Default constructor.
     */
    IBStandardSourceGen();

    /*!
     * \brief Destructor.
     */
    ~IBStandardSourceGen();

    /*!
     * \brief Returns a boolean indicating whether the class has been registered
     * with the singleton IBTK::StreamableManager object.
     */
    static bool getIsRegisteredWithStreamableManager();

    /*!
     * \brief Set the number of internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static void setNumSources(int ln, unsigned int num_sources);

    /*!
     * \brief Get the number of internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static unsigned int getNumSources(int ln);

    /*!
     * \brief Set the names of the internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static void setSourceNames(int ln, const std::vector<std::string>& names);

    /*!
     * \brief Get the names of the internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static const std::vector<std::string>& getSourceNames(int ln);

    /*!
     * \brief Set the sizes of the internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static void setSourceRadii(int ln, const std::vector<double>& radii);

    /*!
     * \brief Get the sizes of the internal sources and sinks on the specified
     * level of the patch hierarchy.
     */
    static const std::vector<double>& getSourceRadii(int ln);

    /*!
     * \brief Return a reference to the vector of source strengths.
     *
     * \note Users \em must \em not change the size of this vector.
     */
    std::vector<double>& getSourceStrengths(int ln);

    /*!
     * \brief Return a const reference to the vector of source strengths.
     */
    const std::vector<double>& getSourceStrengths(int ln) const;

    /*!
     * \brief Return a const reference to the vector of source pressures.
     */
    const std::vector<double>& getSourcePressures(int ln) const;

    /*!
     * \brief Setup the data needed to compute source/sink data on the specified
     * level of the patch hierarchy.
     *
     * \note A default empty implementation is provided.
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
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
    unsigned int getNumSources(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double data_time,
                               IBTK::LDataManager* l_data_manager);

    /*!
     * \brief Compute the source locations for each of the distributed internal
     * sources or sinks.
     *
     * \note Implementations of this method \em must compute the same values for
     * \a X_src on \em each MPI process.  That is to say, \a X_src must provide
     * the location of all of the distributed sources/sinks.
     */
    void getSourceLocations(std::vector<IBTK::Point>& X_src,
                            std::vector<double>& r_src,
                            SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            int level_number,
                            double data_time,
                            IBTK::LDataManager* l_data_manager);

    /*!
     * \brief Set the normalized pressures at the sources.
     */
    void setSourcePressures(const std::vector<double>& P_src,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            int level_number,
                            double data_time,
                            IBTK::LDataManager* l_data_manager);

    /*!
     * \brief Compute the source strengths for each of the distributed internal
     * sources or sinks.
     *
     * \note Implementations of this method \em must compute the same values for
     * \a Q_src on \em each MPI process.  That is to say, \a Q_src must provide
     * the strengths of all of the distributed sources/sinks.
     */
    void computeSourceStrengths(std::vector<double>& Q_src,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                int level_number,
                                double data_time,
                                IBTK::LDataManager* l_data_manager);

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database pointer must be non-null.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBStandardSourceGen(const IBStandardSourceGen& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBStandardSourceGen& operator=(const IBStandardSourceGen& that);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * The numbers of sources/sinks on each level of the patch hierarchy.
     */
    static std::vector<int> s_num_sources;

    /*!
     * The names of the sources and sinks.
     */
    static std::vector<std::vector<std::string> > s_source_names;

    /*!
     * The sizes of the sources and sinks.
     */
    static std::vector<std::vector<double> > s_source_radii;

    /*
     * Source/sink data.
     */
    std::vector<int> d_n_src;
    std::vector<std::vector<std::string> > d_source_names;
    std::vector<std::vector<double> > d_r_src;
    std::vector<std::vector<int> > d_num_perimeter_nodes;
    std::vector<std::vector<double> > d_Q_src, d_P_src;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBStandardSourceGen
