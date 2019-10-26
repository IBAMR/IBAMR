// Filename LSLocateStructureInterface.h
// Created by Nishant Nangia
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
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


/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_LSLocateStructureInterface
#define included_LSLocateStructureInterface

///////////////////////////// INCLUDES ///////////////////////////////////
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/ConstraintIBMethod.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*
 * Pre processing call back function to be hooked into IBAMR:LInitStrategy
 */

void callLSLocateStructureInterfaceCallbackFunction(int D_idx,
                                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                    double time,
                                                    bool initial_time,
                                                    void* ctx);

// Struct to maintain the properties of the rectangle interface
struct RectangleInterface
{
    double length;
    double width;
    double height;
    IBTK::Vector X0;
};

class LSLocateStructureInterface
{
    /*!
     * \brief class LSLocateStructureInterface is a utility class which is used to identify
     * the rectangular interface for level set computations
     */
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateStructureInterface(const std::string& object_name,
                               SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                               SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                               IBTK::LDataManager* lag_data_manager,
                               double vol_elem,
                               RectangleInterface* rectangle);

    /*!
     * Destructor for this class.
     */
    ~LSLocateStructureInterface();

    /*!
     * Reinitialize the level set information
     */
    void setLevelSetPatchData(int D_idx,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              const double time,
                              const bool initial_time);

    //////////////// PRIVATE /////////////////////////////

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    LSLocateStructureInterface();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    LSLocateStructureInterface& operator=(const LSLocateStructureInterface& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    LSLocateStructureInterface(const LSLocateStructureInterface& from);

    /*!
     * Reinitialize the level set information by geometry.
     */
    void setLevelSetPatchDataByGeometry(int D_idx,
                                        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                        const double time,
                                        const bool initial_time);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to the advection-diffusion solver
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variable
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_var;

    /*!
     * IB information
     */
    IBTK::LDataManager* d_lag_data_manager;

    /*!
     * Volume element
     */
    double d_vol_elem;

    /*!
     * Rectangle structure location.
     */
    RectangleInterface* d_rectangle;
};

#endif // #ifndef included_LSLocateStructureInterface
