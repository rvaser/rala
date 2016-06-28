/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <set>

#include "read.hpp"
#include "overlap.hpp"
#include "node.hpp"
#include "edge.hpp"
#include "graph.hpp"

namespace RALAY {

constexpr uint32_t kMaxOverhang = 1000;
constexpr float kMaxOverhangToOverlapRatio = 0.8;

uint32_t classifyOverlap(const std::shared_ptr<Overlap>& ovl) {

    uint32_t overhang = std::min(ovl->a_begin(), ovl->b_begin()) +
        std::min(ovl->a_length() - ovl->a_end(), ovl->b_length() - ovl->b_end());
    uint32_t overlap_length = std::max(ovl->a_end() - ovl->a_begin(),
        ovl->b_end() - ovl->b_begin());

    if (overhang > std::min(kMaxOverhang, (uint32_t) (kMaxOverhangToOverlapRatio * overlap_length))) {
        return 0; // internal match
    }
    if (ovl->a_begin() <= ovl->b_begin() && (ovl->a_length() - ovl->a_end()) <= (ovl->b_length() - ovl->b_end())) {
        return 1; // a contained
    }
    if (ovl->a_begin() >= ovl->b_begin() && (ovl->a_length() - ovl->a_end()) >= (ovl->b_length() - ovl->b_end())) {
        return 2; // b contained
    }
    if (ovl->a_begin() > ovl->b_begin()) {
        return 3; // a to b overlap
    }

    return 4; // b to a overlap
}

std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps) {
    return std::unique_ptr<Graph>(new Graph(reads, overlaps));
}

Graph::Graph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps) :
        nodes_() {

    // remove contained reads and their overlaps before graph construction
    uint32_t max_read_id = 0;
    for (const auto& it: overlaps) {
        max_read_id = std::max(max_read_id, std::max(it->a_id(), it->b_id()));
    }

    std::vector<bool> is_contained(max_read_id + 1, false);
    std::vector<uint8_t> overlap_type(overlaps.size(), 0);

    for (uint32_t i = 0; i < overlaps.size(); ++i) {
        overlap_type[i] = classifyOverlap(overlaps[i]);

        if (overlap_type[i] == 1) { // a contained
            is_contained[overlaps[i]->a_id()] = true;
        } else if (overlap_type[i] == 2) { // b contained
            is_contained[overlaps[i]->b_id()] = true;
        }
    }

    // create assembly graph
    std::vector<int32_t> read_id_to_node_id(max_read_id + 1, -1);
    uint32_t node_id = 0;
    for (const auto& it: reads) {
        if (!is_contained[it->id()]) {
            read_id_to_node_id[it->id()] = node_id;
            nodes_.push_back(createNode(node_id++, it, false)); // normal read
            nodes_.push_back(createNode(node_id++, it, true)); // reverse complement
        }
    }

    uint32_t edge_id = 0;
    for (uint32_t i = 0; i < overlaps.size(); ++i) {
        const auto& it = overlaps[i];
        if (!is_contained[it->a_id()] && !is_contained[it->b_id()]) {
            if (overlap_type[i] < 3) {
                continue;
            }

            uint32_t a = read_id_to_node_id[it->a_id()] + (it->a_rc() == 0 ? 0 : 1);
            uint32_t _a = a + (it->a_rc() == 0 ? 1 : -1);

            uint32_t b = read_id_to_node_id[it->b_id()] + (it->b_rc() == 0 ? 0 : 1);
            uint32_t _b = b + (it->b_rc() == 0 ? 1 : -1);

            if (overlap_type[i] == 3) { // a to b overlap
                edges_.push_back(createEdge(edge_id++, overlaps[i], a, b));
                nodes_[a]->add_suffix_edge_id(edge_id - 1);
                nodes_[b]->add_prefix_edge_id(edge_id - 1);

                edges_.push_back(createEdge(edge_id++, overlaps[i], _b, _a));
                nodes_[_a]->add_prefix_edge_id(edge_id - 1);
                nodes_[_b]->add_suffix_edge_id(edge_id - 1);

            } else if (overlap_type[i] == 4) { // b to a overlap
                edges_.push_back(createEdge(edge_id++, overlaps[i], b, a));
                nodes_[b]->add_suffix_edge_id(edge_id - 1);
                nodes_[a]->add_prefix_edge_id(edge_id - 1);

                edges_.push_back(createEdge(edge_id++, overlaps[i], _a, _b));
                nodes_[_b]->add_prefix_edge_id(edge_id - 1);
                nodes_[_a]->add_suffix_edge_id(edge_id - 1);
            }
        }
    }
}

Graph::~Graph() {
}

void Graph::print() const {

    printf("digraph 1 {\n");
    printf("    overlap = scalexy\n");

    for (const auto& it: nodes_) {

        printf("    %d [label = \"%u: +%u,-%u\"", it->id(), it->id(), it->in_degree(), it->out_degree());
        if (it->id() % 2 == 1) {
            printf(", style = filled, fillcolor = brown1]\n");
            printf("    %d -> %d [style = dotted, arrowhead = none]\n", it->id(), it->id() - 1);
        } else {
            printf("]\n");
        }

        for (const auto& it2: it->suffix_edges_ids()) {
            printf("    %d -> %d", it->id(), edges_[it2]->end_node_id());
            if (it->id() % 2 == 1) {
                printf(" [color = brown1]\n");
            } else {
                printf("\n");
            }
        }
    }

    printf("}\n");
}

}
