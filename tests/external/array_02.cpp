// ---------------------------------------------------------------------
//
// Copyright (c) 2025 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Test that SAMRAI::tbox::Array works as expected. This verifies that the new
// pool-based implementation in IBSAMRAI2 is correct.

#include <ibtk/IBTKInit.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIArray.h>

#include <fstream>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

static std::ostream* output_ptr;
static int ctor_counter = 0;
static int dtor_counter = 0;

// Ensure that we correctly allocate and deallocate nontrivial empty classes
struct A
{
    A()
    {
        (*output_ptr) << "  A::A() " << ctor_counter << "\n";
        ++ctor_counter;
    }

    ~A()
    {
        (*output_ptr) << "  A::~A() " << dtor_counter << "\n";
        ++dtor_counter;
    }
};

template <typename T>
void
test(const std::string& class_name, std::ostream& output)
{
    output << "Test empty tbox::Array<" << class_name << ">\n";
    {
        SAMRAIArray<T> a;
    }

    output << "Test tbox::Array<" << class_name << ">(1)\n";
    {
        SAMRAIArray<T> a(1);
    }

    output << "Test assignment of array tbox::Array<" << class_name << ">\n";
    {
        SAMRAIArray<T> a(1), b(4);
        a = b;
        output << "  Assignment complete\n";
    }

    output << "Test tbox::Array<" << class_name << ">::setNull()\n";
    {
        SAMRAIArray<T> a(4);
        a.setNull();
        output << "  setNull() complete\n";
    }

    output << "Test tbox::Array<" << class_name << ">::resizeArray()\n";
    {
        SAMRAIArray<T> a(4);
        a.resizeArray(3);
        output << "  resizeArray() complete\n";
    }

    // Note that this uses reference counting so we do not see any extra copies being made
    output << "Test tbox::Array<" << class_name << ">::operator=()\n";
    {
        SAMRAIArray<T> a(4), b, c, d;
        b = a;
        c = a;
        d = a;
        TBOX_ASSERT(a.getPointer() == b.getPointer());
        TBOX_ASSERT(a.getPointer() == c.getPointer());
        TBOX_ASSERT(a.getPointer() == d.getPointer());

        output << "  operator=() complete\n";
    }

    // Note that this uses reference counting so we do not see any extra copies being made
    output << "Test tbox::Array<" << class_name << ">::Array(const tbox::Array<" << class_name << "> &)\n";
    {
        SAMRAIArray<T> a(4);
        SAMRAIArray<T> b(a);
        SAMRAIArray<T> c(b);
        SAMRAIArray<T> d(c);

        TBOX_ASSERT(a.getPointer() == b.getPointer());
        TBOX_ASSERT(a.getPointer() == c.getPointer());
        TBOX_ASSERT(a.getPointer() == d.getPointer());

        output << "  'copy' constructor complete\n";
    }
}

int
main()
{
    std::ofstream output("output");
    output_ptr = &output;

    test<A>("A", output);
    test<double>("double", output);
    test<std::string>("std::string", output);
    test<std::vector<double>>("std::vector<double>", output);

    SAMRAIArray<std::string> a(100);
}
