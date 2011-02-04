// Filename: ibtk_utilities.h
// Created on 27 Jan 2011 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_ibtk_utilities
#define included_ibtk_utilities

/////////////////////////////// MACRO DEFINITIONS ////////////////////////////

#include <tbox/PIO.h>

#define IBTK_DO_ONCE(task)                                      \
    do                                                          \
    {                                                           \
        static bool done = false;                               \
        if (done == false)                                      \
        {                                                       \
            task;                                               \
            done = true;                                        \
        }                                                       \
    } while (0)

#define IBTK_DEPRECATED_CLASS1(deprecated_class_name)                   \
    IBTK_DO_ONCE({                                                      \
            SAMRAI::tbox::pout << "WARNING: class "                     \
                               << deprecated_class_name                 \
                               << " is deprecated and may be removed in the future." \
                               << std::endl;                            \
        });

#define IBTK_DEPRECATED_CLASS2(deprecated_class_name,new_class_name)    \
    IBTK_DO_ONCE({                                                      \
            SAMRAI::tbox::pout << "WARNING: class "                     \
                               << deprecated_class_name                 \
                               << " is deprecated and may be removed in the future.\n" \
                               << "Please update your code to use class " \
                               << new_class_name                        \
                               << "."                                   \
                               << std::endl;                            \
        });

#define IBTK_DEPRECATED_FUNCTION1(deprecated_function_name)             \
    IBTK_DO_ONCE({                                                      \
            SAMRAI::tbox::pout << "WARNING: function "                  \
                               << deprecated_function_name              \
                               << " is deprecated and may be removed in the future." \
                               << std::endl;                            \
        });

#define IBTK_DEPRECATED_FUNCTION2(deprecated_function_name,new_function_name) \
    IBTK_DO_ONCE({                                                      \
            SAMRAI::tbox::pout << "WARNING: function "                  \
                               << deprecated_function_name              \
                               << " is deprecated and may be removed in the future.\n" \
                               << "Please update your code to use function " \
                               << new_function_name                     \
                               << "."                                   \
                               << std::endl;                            \
        });

#define IBTK_DEPRECATED_MEMBER_FUNCTION1(class_name,deprecated_method_name) \
    IBTK_DO_ONCE({                                                      \
            SAMRAI::tbox::pout << "WARNING: member function "           \
                               << deprecated_method_name                \
                               << " of class "                          \
                               << class_name                            \
                               << " is deprecated and may be removed in the future." \
                               << std::endl;                            \
        });

#define IBTK_DEPRECATED_MEMBER_FUNCTION2(class_name,deprecated_method_name,new_method_name) \
    IBTK_DO_ONCE({                                                      \
            SAMRAI::tbox::pout << "WARNING: member function "           \
                               << deprecated_method_name                \
                               << " of class "                          \
                               << class_name                            \
                               << " is deprecated and may be removed in the future.\n" \
                               << "Please update your code to use member function " \
                               << new_method_name                       \
                               << "."                                   \
                               << std::endl;                            \
        });

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ibtk_utilities
