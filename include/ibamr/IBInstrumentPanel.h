// Filename: IBInstrumentPanel.h
// Created on 12 May 2007 by Boyce Griffith
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

#ifndef included_IBInstrumentPanel
#define included_IBInstrumentPanel

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <functional>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include "Index.h"
#include "IntVector.h"
#include "boost/multi_array.hpp"
#include "ibtk/ibtk_utilities.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace IBTK
{
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
 * \brief Class IBInstrumentPanel provides support for flow meters and pressure
 * gauges.
 *
 * \note Use of class IBInstrumentPanel requires the Blitz++ array library.
 */
class IBInstrumentPanel : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    IBInstrumentPanel(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~IBInstrumentPanel();

    /*!
     * \return A const reference to the vector of instrument names.
     */
    const std::vector<std::string>& getInstrumentNames() const;

    /*!
     * \return A const reference to the most recent time that the flow meter and
     * pressure gauge values were set.
     *
     * \note This value is not initialized until the first call is made to
     * readInstrumentData().
     */
    const double& getInstrumentDataReadTime() const;

    /*!
     * \return A const reference to the vector of flow meter values.
     *
     * \note This vector is not initialized until the first call is made to
     * readInstrumentData().
     */
    const std::vector<double>& getFlowValues() const;

    /*!
     * \return A const reference to the vector of mean pressure values.
     *
     * \note This vector is not initialized until the first call is made to
     * readInstrumentData().
     */
    const std::vector<double>& getMeanPressureValues() const;

    /*!
     * \return A const reference to the vector of pointwise pressure values.
     *
     * \note This vector is not initialized until the first call is made to
     * readInstrumentData().
     */
    const std::vector<double>& getPointwisePressureValues() const;

    /*!
     * \return A boolean indicating whether there are any instruments embedded
     * within the model data.
     *
     * \note This method returns false and prints a warning message prior to the
     * first call to initializeHierarchyIndependentData().
     */
    bool isInstrumented() const;

    /*!
     * \brief Initialize hierarchy-independent data.
     *
     * The data initialized by this method is assumed \em not to change during
     * the course of a simulation.
     */
    void initializeHierarchyIndependentData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                            IBTK::LDataManager* l_data_manager);

    /*!
     * \brief Initialize hierarchy- and configuration-dependent data.
     */
    void initializeHierarchyDependentData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                          IBTK::LDataManager* l_data_manager,
                                          int timestep_num,
                                          double data_time);

    /*!
     * \brief Compute the flow rates and pressures in the various distributed
     * internal flow meters and pressure gauges.
     */
    void readInstrumentData(int U_data_idx,
                            int P_data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                            IBTK::LDataManager* l_data_manager,
                            int timestep_num,
                            double data_time);

    /*!
     * \brief Set the directory where plot data is to be written.
     *
     * \note By default, visualization data is placed into the directory
     * "viz_inst2d" for two-dimensional simulations and "viz_inst3d" for
     * three-dimensional simulations.
     */
    void setPlotDirectory(const std::string& plot_directory_name);

    /*!
     * \brief Write the plot data to disk.
     */
    void writePlotData(int timestep_num, double simulation_time);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBInstrumentPanel();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBInstrumentPanel(const IBInstrumentPanel& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBInstrumentPanel& operator=(const IBInstrumentPanel& that);

    /*!
     * Read input values, indicated above, from given database.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Output log data to the provided output stream.
     */
    void outputLogData(std::ostream& os);

    /*
     * The object name is used for error reporting purposes.
     */
    std::string d_object_name;

    /*!
     * \brief Instrumentation data.
     */
    bool d_initialized;
    unsigned int d_num_meters;
    std::vector<int> d_num_perimeter_nodes;
    std::vector<IBTK::Vector> d_X_centroid;
    std::vector<boost::multi_array<IBTK::Vector, 1> > d_X_perimeter;
    std::vector<boost::multi_array<IBTK::Vector, 2> > d_X_web, d_dA_web;

    int d_instrument_read_timestep_num;
    double d_instrument_read_time;
    int d_max_instrument_name_len;
    std::vector<std::string> d_instrument_names;
    std::vector<double> d_flow_values, d_mean_pres_values, d_point_pres_values;

    /*!
     * \brief Data structures employed to manage mappings between cell indices
     * and web patch data (i.e., patch centroids and area-weighted normals) and
     * meter centroid data.
     */
    struct IndexFortranOrder : public std::binary_function<SAMRAI::hier::Index<NDIM>, SAMRAI::hier::Index<NDIM>, bool>
    {
        inline bool operator()(const SAMRAI::hier::Index<NDIM>& lhs, const SAMRAI::hier::Index<NDIM>& rhs) const
        {
            return (lhs(0) < rhs(0)
#if (NDIM > 1)
                    ||
                    (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                    ||
                    (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
                        );
        } // operator()
    };

    struct WebPatch
    {
        int meter_num;
        const IBTK::Vector* X;
        const IBTK::Vector* dA;
    };

    typedef std::multimap<SAMRAI::hier::Index<NDIM>, WebPatch, IndexFortranOrder> WebPatchMap;
    std::vector<WebPatchMap> d_web_patch_map;

    struct WebCentroid
    {
        int meter_num;
        const IBTK::Vector* X;
    };

    typedef std::multimap<SAMRAI::hier::Index<NDIM>, WebCentroid, IndexFortranOrder> WebCentroidMap;
    std::vector<WebCentroidMap> d_web_centroid_map;

    /*
     * The directory where data is to be dumped and the most recent timestep
     * number at which data was dumped.
     */
    std::string d_plot_directory_name;

    /*!
     * The log file name and optional flow rate and pressure conversion factors.
     */
    bool d_output_log_file;
    std::string d_log_file_name;
    std::ofstream d_log_file_stream;
    double d_flow_conv, d_pres_conv;
    std::string d_flow_units, d_pres_units;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBInstrumentPanel
