// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/IBTK_MPI.h"
#include "ibtk/StreamableFactory.h"
#include "ibtk/StreamableManager.h"

#include "tbox/Pointer.h"
#include "tbox/ShutdownRegistry.h"

#include <map>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

StreamableManager* StreamableManager::s_data_manager_instance = nullptr;
bool StreamableManager::s_registered_callback = false;
int StreamableManager::s_current_id_number = 0;
const int StreamableManager::s_unregistered_id_number = -1;
unsigned char StreamableManager::s_shutdown_priority = 200;

StreamableManager*
StreamableManager::getManager()
{
    if (!s_data_manager_instance)
    {
        s_data_manager_instance = new StreamableManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_data_manager_instance;
} // getManager

void
StreamableManager::freeManager()
{
    delete s_data_manager_instance;
    s_data_manager_instance = nullptr;
    return;
} // freeManager

int
StreamableManager::getUnregisteredID()
{
    return s_unregistered_id_number;
} // getUnregisteredID

/////////////////////////////// PUBLIC ///////////////////////////////////////

bool
StreamableManager::checkFactoryRegistration(Pointer<StreamableFactory> factory)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(factory);
#endif
    return d_factory_map.count(factory->getStreamableClassID()) == 1;
} // checkFactoryRegistration

int
StreamableManager::registerFactory(Pointer<StreamableFactory> factory)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(factory);
    TBOX_ASSERT(factory->getStreamableClassID() == getUnregisteredID());
#endif
    // These barriers ensure that each factory is assigned the same class ID
    // number on each MPI process.
    IBTK_MPI::barrier();
    const int factory_id = createUniqueID();
    IBTK_MPI::barrier();
    factory->setStreamableClassID(factory_id);
    d_factory_map[factory_id] = factory;
    return factory_id;
} // registerFactory

/////////////////////////////// PROTECTED ////////////////////////////////////

StreamableManager::StreamableManager() : d_factory_map()
{
    // intentionally blank
    return;
} // StreamableManager

StreamableManager::~StreamableManager()
{
    d_factory_map.clear();
    return;
} // ~StreamableManager

int
StreamableManager::createUniqueID()
{
    return s_current_id_number++;
} // createUniqueID

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
