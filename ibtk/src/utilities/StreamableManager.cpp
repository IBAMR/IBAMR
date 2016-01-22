// Filename: StreamableManager.cpp
// Created on 14 Jun 2004 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <map>
#include <ostream>

#include "ibtk/StreamableFactory.h"
#include "ibtk/StreamableManager.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

StreamableManager* StreamableManager::s_data_manager_instance = NULL;
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
    s_data_manager_instance = NULL;
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
    SAMRAI_MPI::barrier();
    const int factory_id = createUniqueID();
    SAMRAI_MPI::barrier();
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
