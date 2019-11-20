// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

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
