/*!
 * @file node.hpp
 *
 * @brief Node class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>

namespace RALAY {

class Read;
class Node;
std::unique_ptr<Node> createNode(uint32_t id, const std::shared_ptr<Read>& read, bool rc);

class Edge;
class Node {
public:

    ~Node();

    uint32_t id() const {
        return id_;
    }

    const std::vector<uint32_t>& prefix_edges_ids() const {
        return prefix_edges_ids_;
    }

    void add_prefix_edge_id(uint32_t id) {
        prefix_edges_ids_.emplace_back(id);
    }

    uint32_t in_degree() const {
        return prefix_edges_ids_.size();
    }

    const std::vector<uint32_t>& suffix_edges_ids() const {
        return suffix_edges_ids_;
    }

    void add_suffix_edge_id(uint32_t id) {
        suffix_edges_ids_.emplace_back(id);
    }

    uint32_t out_degree() const {
        return suffix_edges_ids_.size();
    }

    friend std::unique_ptr<Node> createNode(uint32_t id, const std::shared_ptr<Read>& read, bool rc);

private:

    Node(uint32_t id, const std::shared_ptr<Read>& read, bool rc);
    Node(const Node&) = delete;
    const Node& operator=(const Node&) = delete;

    uint32_t id_;
    std::shared_ptr<Read> read_;
    bool rc_;
    std::vector<uint32_t> prefix_edges_ids_;
    std::vector<uint32_t> suffix_edges_ids_;
};

}
