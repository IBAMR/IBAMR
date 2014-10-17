// Filename: LInitStrategy.h
// Created on 11 Jul 2004 by Boyce Griffith
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

#ifndef included_LInitStrategy
#define included_LInitStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <string>
#include <utility>

#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

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
} // namespace SAMRAI

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class LDataManager;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LInitStrategy provides a mechanism for specifying the
 * initial configuration of the curvilinear mesh.
 */
class LInitStrategy : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    LInitStrategy();

    /*!
     * \brief Destructor.
     */
    virtual ~LInitStrategy();

    /*!
     * \return A boolean value indicating whether Lagrangian data is associated
     * with the given level in the patch hierarchy.
     */
    virtual bool getLevelHasLagrangianData(int level_number, bool can_be_refined) const = 0;

    /*!
     * \return The number of local nodes on the patch level.
     */
    virtual unsigned int
    computeLocalNodeCountOnPatchLevel(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
                                      double init_data_time,
                                      bool can_be_refined,
                                      bool initial_time) = 0;

    /*!
     * \brief Initialize the structure indexing information on the patch level.
     *
     * \note A default empty implementation is provided.
     */
    virtual void
    initializeStructureIndexingOnPatchLevel(std::map<int, std::string>& strct_id_to_strct_name_map,
                                            std::map<int, std::pair<int, int> >& strct_id_to_lag_idx_range_map,
                                            int level_number,
                                            double init_data_time,
                                            bool can_be_refined,
                                            bool initial_time,
                                            LDataManager* l_data_manager);

    /*!
     * \brief Initialize the LNode and LData data needed to specify the
     * configuration of the curvilinear mesh on the patch level.
     *
     * \return The number of local nodes initialized on the patch level.
     */
    virtual unsigned int
    initializeDataOnPatchLevel(int lag_node_index_idx,
                               unsigned int global_index_offset,
                               unsigned int local_index_offset,
                               SAMRAI::tbox::Pointer<LData> X_data,
                               SAMRAI::tbox::Pointer<LData> U_data,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double init_data_time,
                               bool can_be_refined,
                               bool initial_time,
                               LDataManager* l_data_manager) = 0;

    /*!
     * \brief Initialize the LData needed to specify the mass and spring
     * constant data required by the penalty IB method.
     *
     * \return The number of local nodes initialized on the patch level.
     *
     * \note A default empty implementation is provided when support for massive
     * boundaries is not required.
     */
    virtual unsigned int
    initializeMassDataOnPatchLevel(unsigned int global_index_offset,
                                   unsigned int local_index_offset,
                                   SAMRAI::tbox::Pointer<LData> M_data,
                                   SAMRAI::tbox::Pointer<LData> K_data,
                                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                   int level_number,
                                   double init_data_time,
                                   bool can_be_refined,
                                   bool initial_time,
                                   LDataManager* l_data_manager);

    /*!
     * \brief Initialize the LData needed to specify director vectors required
     * by some material models.
     *
     * \return The number of local nodes initialized on the patch level.
     *
     * \note A default empty implementation is provided when support for
     * directors is not required.
     */
    virtual unsigned int
    initializeDirectorDataOnPatchLevel(unsigned int global_index_offset,
                                       unsigned int local_index_offset,
                                       SAMRAI::tbox::Pointer<LData> D_data,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       int level_number,
                                       double init_data_time,
                                       bool can_be_refined,
                                       bool initial_time,
                                       LDataManager* l_data_manager);

    /*!
     * \brief Provide cell tagging for the initial configuration of the
     * Lagrangian mesh.
     *
     * When the patch hierarchy is being constructed at the initial simulation
     * time, it is necessary that the gridding algorithm be instructed where to
     * place local refinement in order to accommodate portions of the curvilinear
     * mesh that will reside in the yet-to-be-constructed level(s) of the patch
     * hierarchy.
     *
     * \note A default empty implementation is provided when support for local
     * mesh refinement is not required.
     */
    virtual void tagCellsForInitialRefinement(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                              int level_number,
                                              double error_data_time,
                                              int tag_index);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LInitStrategy(const LInitStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LInitStrategy& operator=(const LInitStrategy& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LInitStrategy
