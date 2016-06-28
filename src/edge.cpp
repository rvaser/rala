/*!
 * @file edge.hpp
 *
 * @brief Edge class source file
 */

#include "overlap.hpp"
#include "edge.hpp"

namespace RALAY {

std::unique_ptr<Edge> createEdge(uint32_t id, const std::shared_ptr<Overlap>& overlap,
    uint32_t begin_node_id, uint32_t end_node_id) {

    return std::unique_ptr<Edge>(new Edge(id, overlap, begin_node_id, end_node_id));
}

Edge::Edge(uint32_t id, const std::shared_ptr<Overlap>& overlap, uint32_t begin_node_id,
    uint32_t end_node_id) :
        id_(id), overlap_(overlap), begin_node_id_(begin_node_id), end_node_id_(end_node_id) {
}

Edge::~Edge() {
}

}
