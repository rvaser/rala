/*!
 * @file node.cpp
 *
 * @brief Node class source file
 */

#include "read.hpp"
#include "node.hpp"

namespace RALAY {

std::unique_ptr<Node> createNode(uint32_t id, const std::shared_ptr<Read>& read, bool rc) {
    return std::unique_ptr<Node>(new Node(id, read, rc));
}

Node::Node(uint32_t id, const std::shared_ptr<Read>& read, bool rc) :
        id_(id), read_(read), rc_(rc), prefix_edges_ids_(), suffix_edges_ids_() {
}

Node::~Node() {
}

}
