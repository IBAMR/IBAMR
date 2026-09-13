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

#ifndef included_VelocityInitialCondition
#define included_VelocityInitialCondition

#include <ibtk/CartGridFunction.h>

#include <PointwiseLevelSet.h>

#include <string>
#include <vector>

class VelocityInitialCondition : public IBTK::CartGridFunction
{
public:
    VelocityInitialCondition(const std::string& object_name,
                             double num_interface_cells,
                             std::vector<double> inside_velocity,
                             std::vector<double> outside_velocity,
                             SAMRAI::tbox::Pointer<MultiphaseExamples::SphereLevelSet> sphere);
    bool isTimeDependent() const override;
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level = nullptr) override;

private:
    double d_num_interface_cells;
    std::vector<double> d_inside_velocity, d_outside_velocity;
    SAMRAI::tbox::Pointer<MultiphaseExamples::SphereLevelSet> d_sphere;
};

#endif
