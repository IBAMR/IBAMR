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

#include <CoupledApplication.h>

#include <ibamr/app_namespaces.h>

namespace
{
class StefanExpansion : public PhaseChangeExamples::CoupledApplication
{
public:
    explicit StefanExpansion(Pointer<AppInitializer> app_initializer) : CoupledApplication(app_initializer)
    {
        useEnthalpy();
        registerLevelSet();
        registerLiquidFraction();
        registerLiquidFractionGradient();
        registerSpecificEnthalpy();
        registerHeaviside();
        registerTransportFields();
        registerMaterialFields();
        registerLiquidFractionTagging();
        setHeavisideBoundary();
        setTemperatureBoundary();
        setEnthalpyBoundary();
        setLiquidFractionBoundary();
        setConductivityBoundary();
        setVelocityBoundary();
        setDensityBoundary();
        setViscosityBoundary();
        setLevelSetBoundary();
        registerMaterialProperties();
        registerPhaseChangeSources("F");
        registerSolidDrag();
    }
};

std::unique_ptr<PhaseChangeExamples::CoupledApplication>
create_application(Pointer<AppInitializer> app_initializer)
{
    return std::make_unique<StefanExpansion>(app_initializer);
}
} // namespace

int
main(int argc, char* argv[])
{
    return PhaseChangeExamples::run_coupled(argc, argv, &create_application);
}
