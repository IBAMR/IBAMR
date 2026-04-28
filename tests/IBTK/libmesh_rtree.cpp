#include <ibtk/libmesh_rtree_wrappers.h>

#include <libmesh/elem.h>
#include <libmesh/libmesh.h>
#include <libmesh/libmesh_config.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/replicated_mesh.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <fstream>

// verify that the rtree wrappers work as intended.

int
main(int argc, char** argv)
{
    libMesh::LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    libMesh::ReplicatedMesh mesh(init.comm(), 3);
    libMesh::MeshTools::Generation::build_sphere(mesh, 1.0, 1, libMesh::HEX8);

    std::ofstream output("output");

    std::vector<libMesh::BoundingBox> bounding_boxes;
    std::vector<libMesh::Elem*> elems;
    const auto el_begin = mesh.local_elements_begin();
    const auto el_end = mesh.local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        elems.emplace_back(*el_it);
        bounding_boxes.emplace_back((*el_it)->loose_bounding_box());
    }

    const auto rtree = IBTK::pack_rtree_of_indices(bounding_boxes);

    output << "Number of elements = " << mesh.n_elem() << std::endl;

    // loop over element bounding boxes...
    for (unsigned int i = 0; i < bounding_boxes.size(); ++i)
    {
        // and determine which patches each intersects.
        auto print_neighbor = [&](const std::size_t& elem_n) {
            output << "  near element " << elems[elem_n]->id() << " with centroid " << elems[elem_n]->true_centroid()
                   << '\n';
        };

        output << "element " << elems[i]->id() << " with centroid " << elems[i]->true_centroid() << '\n';
        namespace bgi = boost::geometry::index;
        rtree.query(bgi::intersects(bounding_boxes[i]), boost::make_function_output_iterator(print_neighbor));
    }
}
