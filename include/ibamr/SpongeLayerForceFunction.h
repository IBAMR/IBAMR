// Filename: SpongeLayerForceFunction.h
// Created on 28 Oct 2011 by Boyce Griffith
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

#ifndef included_SpongeLayerForceFunction
#define included_SpongeLayerForceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "boost/array.hpp"
#include "ibtk/CartGridFunction.h"
#include "SAMRAI/tbox/Array.h"

namespace IBAMR
{
class INSHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{

class Variable;

class Patch;
} // namespace hier
namespace pdat
{
template <class TYPE>
class SideData;
template <class TYPE>
class CellData;
} // namespace pdat
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class SpongeLayerForceFunction provides forcing at physical boundaries
 * that weakly imposes homogeneous Dirichlet boundary conditions.
 */
class SpongeLayerForceFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    SpongeLayerForceFunction(const std::string& object_name,
                             boost::shared_ptr<SAMRAI::tbox::Database> input_db,
                             const INSHierarchyIntegrator* fluid_solver,
                             boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry);

    /*!
     * \brief Destructor.
     */
    ~SpongeLayerForceFunction();

    /*!
     * \name Methods to set the data.
     */
    //\{

    /*!
     * \note This concrete IBTK::CartGridFunction is time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        boost::shared_ptr<SAMRAI::hier::Variable> var,
                        boost::shared_ptr<SAMRAI::hier::Patch> patch,
                        double data_time,
                        bool initial_time = false,
                        boost::shared_ptr<SAMRAI::hier::PatchLevel> level = NULL);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SpongeLayerForceFunction();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SpongeLayerForceFunction(const SpongeLayerForceFunction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SpongeLayerForceFunction& operator=(const SpongeLayerForceFunction& that);

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatchCell(boost::shared_ptr<SAMRAI::pdat::CellData<double> > F_data,
                            boost::shared_ptr<SAMRAI::pdat::CellData<double> > U_current_data,
                            boost::shared_ptr<SAMRAI::pdat::CellData<double> > U_new_data,
                            double kappa,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch);

    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatchSide(boost::shared_ptr<SAMRAI::pdat::SideData<double> > F_data,
                            boost::shared_ptr<SAMRAI::pdat::SideData<double> > U_current_data,
                            boost::shared_ptr<SAMRAI::pdat::SideData<double> > U_new_data,
                            double kappa,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch);

    boost::array<std::vector<bool>, 2 * NDIM> d_forcing_enabled;
    boost::array<double, 2 * NDIM> d_width;
    const INSHierarchyIntegrator* const d_fluid_solver;
    boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> d_grid_geometry;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SpongeLayerForceFunction
