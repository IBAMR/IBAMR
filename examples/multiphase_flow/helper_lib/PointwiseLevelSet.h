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

#ifndef included_MultiphaseExamples_PointwiseLevelSet
#define included_MultiphaseExamples_PointwiseLevelSet

#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

#include <string>

namespace MultiphaseExamples
{
// Analytic geometry evaluation is separate from LSLocateInterface's initialize-then-copy policy.
// Reuse CartGridFunction's hierarchy and level traversal for every pointwise geometry.
class PointwiseLevelSet : public IBTK::CartGridFunction
{
public:
    explicit PointwiseLevelSet(const std::string& object_name);
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level = nullptr) override;
    virtual double evaluateSignedDistance(const IBTK::Vector& X, double time) const = 0;
};

class SphereLevelSet : public PointwiseLevelSet
{
public:
    SphereLevelSet(const std::string& object_name, const IBTK::Vector& center, double radius);
    bool isTimeDependent() const override;
    double evaluateSignedDistance(const IBTK::Vector& X, double time) const override;

private:
    IBTK::Vector d_center;
    double d_radius;
};

// Coordinate-minus-height orientation, as used by the multiphase layer example.
class PlaneLevelSet : public PointwiseLevelSet
{
public:
    PlaneLevelSet(const std::string& object_name, int axis, double height);
    bool isTimeDependent() const override;
    double evaluateSignedDistance(const IBTK::Vector& X, double time) const override;

private:
    int d_axis;
    double d_height;
};
} // namespace MultiphaseExamples

#endif
