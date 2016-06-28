/*!
 * @file edge.hpp
 *
 * @brief Edge class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>

namespace RALAY {

class Overlap;
class Node;
class Edge;
std::unique_ptr<Edge> createEdge(uint32_t id, const std::shared_ptr<Overlap>& overlap,
    uint32_t begin_node_id, uint32_t end_node_id);

class Edge {
public:

    ~Edge();

    uint32_t id() const {
        return id_;
    }

    uint32_t begin_node_id() const {
        return begin_node_id_;
    }

    uint32_t end_node_id() const {
        return end_node_id_;
    }

    friend std::unique_ptr<Edge> createEdge(uint32_t id, const std::shared_ptr<Overlap>& overlap,
        uint32_t begin_node_id, uint32_t end_node_id);

private:

    Edge(uint32_t id, const std::shared_ptr<Overlap>& overlap, uint32_t begin_node_id,
        uint32_t end_node_id);
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    uint32_t id_;
    std::shared_ptr<Overlap> overlap_;
    uint32_t begin_node_id_;
    uint32_t end_node_id_;
};

}
