// Filename: SimpleIBSourceGen.C
// Created on 01 Jun 2009 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "SimpleIBSourceGen.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

SimpleIBSourceGen::SimpleIBSourceGen()
{
    if (NDIM != 2)
    {
        TBOX_ERROR("This implementation is for NDIM == 2 only.\n");
    }
    return;
}// SimpleIBSourceGen

SimpleIBSourceGen::~SimpleIBSourceGen()
{
    // intentionally blank
    return;
}// ~SimpleIBSourceGen

void
SimpleIBSourceGen::setTimeInterval(
    const double current_time,
    const double new_time)
{
    // Method setTimeInterval() is used to report to class SimpleIBSourceGen the
    // current simulation time.  This simple implementation of the strategy
    // class IBLagrangianSourceStrategy does not make use of this feature,
    // however.
    return;
}// setTimeInterval

void
SimpleIBSourceGen::initializeLevelData(
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    // Method initializeLevelData() is used to initialize any hierarchy-specific
    // data required to compute the source/sink strengths and is called after
    // any adaptive regridding operations.  This simple implementation of the
    // strategy class IBLagrangianSourceStrategy does not make use of this
    // feature, however.
    return;
}// initializeLevelData

int
SimpleIBSourceGen::getNumSources(
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    // Method getNumSources() is used to specify the number of sources/sinks on
    // each level of the patch hierarchy.  In general, you probably want
    // sources/sinks on the finest level; however, in this simple
    // implementation, we simply indicate that there are four sources on the
    // coarsest level of the patch hierarchy.
    //
    // NOTE: Class IBHierarchyIntegrator is setup to tag cells for refinement in
    // the vicinity of fluid sources/sinks.  This ensures that sources/sinks are
    // surrounded by the finest level of the AMR grid.
    if (level_number == 0)
    {
        return 4;
    }
    else
    {
        return 0;
    }
}// getNumSources

void
SimpleIBSourceGen::getSourceLocations(
    std::vector<std::vector<double> >& X_src,
    std::vector<double>& r_src,
    tbox::Pointer<LNodeLevelData> X_data,
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    // Recall that we have 4 sources/sinks located on level 0 of the patch
    // hierarchy.
    if (level_number != 0) return;

    // Set the positions of the 4 sources on level 0.
    //
    // NOTE: Because X_data is provided as an argument to this function, it is
    // possible to set the positions of the sources in terms of the
    // configuration of the immersed elastic structure.
    TBOX_ASSERT(X_src.size() == 4);
    X_src[0][0] = 0.25; X_src[0][1] = 0.25;  // src 0 is located at (0.25,0.25)
    X_src[1][0] = 0.75; X_src[1][1] = 0.25;  // src 1 is located at (0.75,0.25)
    X_src[2][0] = 0.25; X_src[2][1] = 0.75;  // src 2 is located at (0.25,0.75)
    X_src[3][0] = 0.75; X_src[3][1] = 0.75;  // src 3 is located at (0.75,0.75)

    // Set the radii of the 4 sources on level 0.
    //
    // NOTE: Each source radius is automatically reset in IBHierarchyIntegrator
    // to be the closest integer multiple of the grid spacing which is larger
    // than the prescribed radius.  This is required by the shape function used
    // to "spread out" the point sources/sinks.
    TBOX_ASSERT(r_src.size() == 4);
    r_src[0] = 4.0/64.0;
    r_src[1] = 8.0/64.0;
    r_src[2] = 8.0/64.0;
    r_src[3] = 4.0/64.0;
    return;
}// getSourceLocations

void
SimpleIBSourceGen::setSourcePressures(
    const std::vector<double>& P_src,
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    // Method setSourcePressures() is used to provide to class SimpleIBSourceGen
    // the pressures at each of the internal sources/sinks relative to the mean
    // pressure in the external source sink (see below).  The pressure data may
    // be used to compute pressure-dependent inflow rates in method
    // computeSourceStrengths(), but this simple implementation of the strategy
    // class IBLagrangianSourceStrategy does not make use of this feature.
    return;
}// setSourcePressures

void
SimpleIBSourceGen::computeSourceStrengths(
    std::vector<double>& Q_src,
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    // Recall that we have 4 sources/sinks located on level 0 of the patch
    // hierarchy.
    if (level_number != 0) return;

    // Set the source strengths of the 4 sources on level 0.
    //
    // NOTE: The net flow through the sources and sinks DOES NOT need to be
    // zero.  An external source/sink is automatically generated by class
    // IBHierarchyIntegrator along the boundary of the computational domain to
    // ensure that there is no net inflow into the fluid volume (a requirement
    // to obtain discrete solvability!).  This external source/sink is always
    // placed along the "upper" and "lower" boundaries, e.g., in 2D, the
    // top/bottom y boundaries, and in 3D, the top/bottom z boundaries.
    TBOX_ASSERT(Q_src.size() == 4);
    Q_src[0] = +1.0;
    Q_src[1] = -1.0;
    Q_src[2] = +1.0;
    Q_src[3] = -1.0;
    return;
}// computeSourceStrengths

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
