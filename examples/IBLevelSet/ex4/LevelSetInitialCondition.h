// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_LevelSetInitialCondition
#define included_LevelSetInitialCondition

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <ibamr/app_namespaces.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

// Struct to maintain the properties of the circular interface
struct FoilInterface
{
    Eigen::Vector3d X0, X1, X2, X3;
    double R, mass;
    double freq, theta_0;
};

/*!
 * \brief Class LevelSetInitialCondition provides an initial condition for
 * the level set function.
 */
class LevelSetInitialCondition : public CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    LevelSetInitialCondition(const std::string& object_name, FoilInterface init_foil);

    /*!
     * \brief Empty destructor.
     */
    ~LevelSetInitialCondition();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete LevelSetInitialCondition object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        Pointer<Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL));

    //\}

    inline int sgn(double v)
    {
        return ((v < 0) ? -1 : (v > 0) ? 1 : 0);
    }

private:
    LevelSetInitialCondition();

    LevelSetInitialCondition(const LevelSetInitialCondition& from);

    LevelSetInitialCondition& operator=(const LevelSetInitialCondition& that);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Initial level set information.
     */
    FoilInterface d_init_foil;
};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LevelSetInitialCondition
