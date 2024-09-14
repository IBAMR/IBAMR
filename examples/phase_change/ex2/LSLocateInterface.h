// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
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

#ifndef included_LSLocateInterface
#define included_LSLocateInterface

///////////////////////////// INCLUDES ///////////////////////////////////
#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*
 * Pre processing call back function to be hooked into IBAMR:LInitStrategy
 */

void callLSLocateInterfaceCallbackFunction(int D_idx,
                                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                           double time,
                                           bool initial_time,
                                           void* ctx);

class LSLocateInterface
{
    /*!
     * \brief class LSLocateInterface is a utility class which is used to identify
     * the nterface for level set computations
     */
public:
    /*!
     * The only constructor of this class.
     */
    LSLocateInterface(const std::string& object_name,
                      SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                      const double bubble_radius);

    /*!
     * Destructor for this class.
     */
    ~LSLocateInterface();

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
    LSLocateInterface();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    LSLocateInterface& operator=(const LSLocateInterface& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    LSLocateInterface(const LSLocateInterface& from);

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
     * Radius of a bubble
     */
    double d_radius;
};

#endif // #ifndef included_LSLocateInterface
