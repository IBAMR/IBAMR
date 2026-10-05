// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Check that a tbox::ParallelBuffer, the stream buffer of tbox::pout,
// tbox::perr, and tbox::plog, forwards the characters that it is given,
// including null characters, and reads no others.

#include <sys/mman.h>
#include <tbox/ParallelBuffer.h>
#include <tbox/Utilities.h>

#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>

#include <ibtk/app_namespaces.h>

int
main()
{
    std::ofstream output("output");

    // Characters are written from the end of a page of memory that is followed
    // by an inaccessible page, so that reading past the last character fails.
    const std::size_t page_size = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
    void* const pages = mmap(nullptr, 2 * page_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    TBOX_ASSERT(pages != MAP_FAILED);
    char* const page_end = static_cast<char*>(pages) + page_size;
    const int ierr = mprotect(page_end, page_size, PROT_NONE);
    TBOX_ASSERT(ierr == 0);

    std::ostringstream received;
    tbox::ParallelBuffer buffer;
    buffer.setOutputStream1(&received);
    std::ostream stream(&buffer);
    std::string written;
    const auto write = [&](const std::string& chars)
    {
        std::copy(chars.begin(), chars.end(), page_end - chars.size());
        stream.write(page_end - chars.size(), static_cast<std::streamsize>(chars.size()));
        written += chars;
    };

    // A line, and then a line that contains a null character and is long
    // enough that the buffer grows while it holds the null character.
    write("ab\n");
    write(std::string("c\0d", 3));
    write(std::string(1000, 'e') + "\n");

    const std::string result = received.str();
    std::size_t n_incorrect = 0;
    for (std::size_t k = 0; k < std::min(written.size(), result.size()); ++k)
    {
        if (result[k] != written[k])
        {
            ++n_incorrect;
        }
    }
    output << "characters written: " << written.size() << '\n';
    output << "characters received: " << result.size() << '\n';
    output << "incorrect characters: " << n_incorrect << '\n';

    munmap(pages, 2 * page_size);
}
