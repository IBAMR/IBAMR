// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/IBTKInit.h>
#include <ibtk/box_utilities.h>

#include <fstream>
#include <vector>

using namespace SAMRAI;
// check that @p box is completely contained by one of the boxes in @p boxes.
bool
check_box_contained(const std::vector<hier::Box<NDIM> >& boxes, const hier::Box<NDIM>& box)
{
    for (const auto& b : boxes)
        if (b.contains(box)) return true;
    return false;
}

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    std::ofstream out("output");

    if (NDIM == 2)
    {
        {
            out << "Box merge result 1\n";
            std::vector<hier::Box<NDIM> > boxes;
            boxes.emplace_back(hier::Index<NDIM>(192, 140), hier::Index<NDIM>(263, 199));
            boxes.emplace_back(hier::Index<NDIM>(192, 200), hier::Index<NDIM>(239, 223));
            boxes.emplace_back(hier::Index<NDIM>(192, 224), hier::Index<NDIM>(215, 255));
            boxes.emplace_back(hier::Index<NDIM>(240, 200), hier::Index<NDIM>(263, 223));
            boxes.emplace_back(hier::Index<NDIM>(216, 224), hier::Index<NDIM>(239, 255));
            boxes.emplace_back(hier::Index<NDIM>(240, 224), hier::Index<NDIM>(263, 255));

            const auto result = IBTK::merge_boxes_by_longest_edge(boxes);
            for (const auto& box : result) out << box << '\n';
            for (const auto& box : boxes) check_box_contained(result, box);
        }

        {
            out << "Box merge result 2\n";
            std::vector<hier::Box<NDIM> > boxes;
            boxes.emplace_back(hier::Index<NDIM>(264, 140), hier::Index<NDIM>(307, 179));
            boxes.emplace_back(hier::Index<NDIM>(264, 180), hier::Index<NDIM>(307, 203));
            boxes.emplace_back(hier::Index<NDIM>(264, 204), hier::Index<NDIM>(307, 227));
            boxes.emplace_back(hier::Index<NDIM>(264, 228), hier::Index<NDIM>(307, 255));

            const auto result = IBTK::merge_boxes_by_longest_edge(boxes);
            for (const auto& box : result) out << box << '\n';
            for (const auto& box : boxes) check_box_contained(result, box);
        }

        {
            out << "Box merge result 3\n";
            std::vector<hier::Box<NDIM> > boxes;
            boxes.emplace_back(hier::Index<NDIM>(268, 256), hier::Index<NDIM>(307, 279));
            boxes.emplace_back(hier::Index<NDIM>(192, 256), hier::Index<NDIM>(215, 287));
            boxes.emplace_back(hier::Index<NDIM>(192, 288), hier::Index<NDIM>(243, 311));
            boxes.emplace_back(hier::Index<NDIM>(216, 256), hier::Index<NDIM>(239, 287));
            boxes.emplace_back(hier::Index<NDIM>(240, 256), hier::Index<NDIM>(267, 287));
            boxes.emplace_back(hier::Index<NDIM>(244, 288), hier::Index<NDIM>(267, 311));

            const auto result = IBTK::merge_boxes_by_longest_edge(boxes);
            for (const auto& box : result) out << box << '\n';
            for (const auto& box : boxes) check_box_contained(result, box);
        }

        {
            out << "Box merge result 4\n";
            std::vector<hier::Box<NDIM> > boxes;
            boxes.emplace_back(hier::Index<NDIM>(348, 292), hier::Index<NDIM>(391, 315));
            boxes.emplace_back(hier::Index<NDIM>(308, 332), hier::Index<NDIM>(347, 371));
            boxes.emplace_back(hier::Index<NDIM>(376, 256), hier::Index<NDIM>(399, 291));
            boxes.emplace_back(hier::Index<NDIM>(400, 256), hier::Index<NDIM>(423, 291));
            boxes.emplace_back(hier::Index<NDIM>(392, 292), hier::Index<NDIM>(423, 315));
            boxes.emplace_back(hier::Index<NDIM>(348, 316), hier::Index<NDIM>(423, 371));

            const auto result = IBTK::merge_boxes_by_longest_edge(boxes);
            for (const auto& box : result) out << box << '\n';
            for (const auto& box : boxes) check_box_contained(result, box);
        }

        {
            // We should successfully do nothing with an empty vector of boxes
            const auto result_2 = IBTK::merge_boxes_by_longest_edge({});
            TBOX_ASSERT(result_2.size() == 0);
        }
    }
    if (NDIM == 3)
    {
        {
            out << "Box merge result 1\n";
            std::vector<hier::Box<NDIM> > boxes;
            boxes.emplace_back(hier::Index<NDIM>(36, 24, 24), hier::Index<NDIM>(59, 63, 63));
            boxes.emplace_back(hier::Index<NDIM>(60, 24, 24), hier::Index<NDIM>(75, 47, 63));
            boxes.emplace_back(hier::Index<NDIM>(60, 48, 24), hier::Index<NDIM>(75, 63, 47));
            boxes.emplace_back(hier::Index<NDIM>(60, 48, 48), hier::Index<NDIM>(75, 63, 63));
            boxes.emplace_back(hier::Index<NDIM>(36, 24, 64), hier::Index<NDIM>(59, 63, 103));
            boxes.emplace_back(hier::Index<NDIM>(60, 24, 64), hier::Index<NDIM>(75, 47, 103));
            boxes.emplace_back(hier::Index<NDIM>(60, 48, 64), hier::Index<NDIM>(75, 63, 79));
            boxes.emplace_back(hier::Index<NDIM>(36, 64, 24), hier::Index<NDIM>(59, 103, 63));
            boxes.emplace_back(hier::Index<NDIM>(60, 64, 24), hier::Index<NDIM>(75, 79, 47));
            boxes.emplace_back(hier::Index<NDIM>(60, 64, 48), hier::Index<NDIM>(75, 79, 63));
            boxes.emplace_back(hier::Index<NDIM>(36, 64, 64), hier::Index<NDIM>(59, 103, 103));
            boxes.emplace_back(hier::Index<NDIM>(60, 64, 64), hier::Index<NDIM>(75, 79, 79));
            boxes.emplace_back(hier::Index<NDIM>(76, 24, 24), hier::Index<NDIM>(91, 47, 63));
            boxes.emplace_back(hier::Index<NDIM>(76, 48, 24), hier::Index<NDIM>(91, 63, 47));
            boxes.emplace_back(hier::Index<NDIM>(76, 48, 48), hier::Index<NDIM>(91, 63, 63));
            boxes.emplace_back(hier::Index<NDIM>(76, 24, 64), hier::Index<NDIM>(91, 47, 103));
            boxes.emplace_back(hier::Index<NDIM>(76, 48, 64), hier::Index<NDIM>(91, 63, 79));
            boxes.emplace_back(hier::Index<NDIM>(76, 64, 24), hier::Index<NDIM>(91, 79, 47));
            boxes.emplace_back(hier::Index<NDIM>(76, 64, 48), hier::Index<NDIM>(91, 79, 63));
            boxes.emplace_back(hier::Index<NDIM>(76, 64, 64), hier::Index<NDIM>(91, 79, 79));
            boxes.emplace_back(hier::Index<NDIM>(60, 48, 80), hier::Index<NDIM>(75, 63, 103));
            boxes.emplace_back(hier::Index<NDIM>(60, 64, 80), hier::Index<NDIM>(75, 79, 103));
            boxes.emplace_back(hier::Index<NDIM>(76, 48, 80), hier::Index<NDIM>(91, 63, 103));
            boxes.emplace_back(hier::Index<NDIM>(76, 64, 80), hier::Index<NDIM>(91, 79, 103));
            boxes.emplace_back(hier::Index<NDIM>(60, 80, 24), hier::Index<NDIM>(75, 103, 63));
            boxes.emplace_back(hier::Index<NDIM>(60, 80, 64), hier::Index<NDIM>(75, 103, 103));
            boxes.emplace_back(hier::Index<NDIM>(76, 80, 24), hier::Index<NDIM>(91, 103, 63));
            boxes.emplace_back(hier::Index<NDIM>(76, 80, 64), hier::Index<NDIM>(91, 103, 103));
            boxes.emplace_back(hier::Index<NDIM>(92, 24, 24), hier::Index<NDIM>(115, 63, 63));
            boxes.emplace_back(hier::Index<NDIM>(92, 24, 64), hier::Index<NDIM>(115, 63, 103));
            boxes.emplace_back(hier::Index<NDIM>(92, 64, 24), hier::Index<NDIM>(115, 103, 63));
            boxes.emplace_back(hier::Index<NDIM>(92, 64, 64), hier::Index<NDIM>(115, 103, 103));
            const auto result = IBTK::merge_boxes_by_longest_edge(boxes);
            for (const auto& box : result) out << box << '\n';
            for (const auto& box : boxes)
            {
                TBOX_ASSERT(check_box_contained(result, box));
            }
        }

        {
            // We should successfully do nothing with an empty vector of boxes
            const auto result_2 = IBTK::merge_boxes_by_longest_edge({});
            TBOX_ASSERT(result_2.size() == 0);
        }
    }
}
