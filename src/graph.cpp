/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <list>

#include "read.hpp"
#include "overlap.hpp"
#include "graph.hpp"

namespace RALAY {

constexpr uint32_t kMaxOverhang = 1000;
constexpr float kMaxOverhangToOverlapRatio = 0.8;

void trimReads(std::vector<std::shared_ptr<Read>>& reads, std::vector<std::shared_ptr<Overlap>>& overlaps) {

}

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

class Graph::Node {
public:

    Node(uint32_t _id, const std::shared_ptr<Read>& _read) :
        id(_id), pair(), read(_read), prefix_edges(), suffix_edges(),
        rc(id % 2 == 1), mark(false) {
    }
    Node(const Node&) = delete;
    const Node& operator=(const Node&) = delete;

    ~Node() {}

    uint32_t read_id() const {
        return read->id();
    }

    uint32_t length() const {
        return read->sequence().size();
    }

    uint32_t in_degree() const {
        return prefix_edges.size();
    }

    uint32_t out_degree() const {
        return suffix_edges.size();
    }

    uint32_t id;
    Node* pair;
    std::shared_ptr<Read> read;
    std::list<Edge*> prefix_edges;
    std::list<Edge*> suffix_edges;
    bool rc;
    bool mark;
};

class Graph::Edge {
public:

    Edge(uint32_t _id, const std::shared_ptr<Overlap>& _overlap, Node* _begin_node,
        Node* _end_node) :
            id(_id), pair(), overlap(_overlap), begin_node(_begin_node),
            end_node(_end_node), mark(false) {
    }
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    ~Edge() {}

    uint32_t length() const {
        if (begin_node->read_id() == overlap->a_id()) {
            return overlap->a_begin();
        }
        return overlap->b_begin();
    }

    uint32_t id;
    Edge* pair;
    std::shared_ptr<Overlap> overlap;
    Node* begin_node;
    Node* end_node;
    bool mark;
};

std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps) {
    return std::unique_ptr<Graph>(new Graph(reads, overlaps));
}

Graph::Graph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps) :
        nodes_(), edges_() {

    // remove contained reads and their overlaps before graph construction
    uint32_t max_read_id = 0;
    for (const auto& it: reads) {
        max_read_id = std::max(max_read_id, it->id());
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

            Node* node = new Node(node_id++, it); // normal read
            Node* _node = new Node(node_id++, it); // reverse complement

            node->pair = _node;
            _node->pair = node;

            nodes_.push_back(std::unique_ptr<Node>(node));
            nodes_.push_back(std::unique_ptr<Node>(_node));
        }
    }

    uint32_t edge_id = 0;
    for (uint32_t i = 0; i < overlaps.size(); ++i) {
        const auto& it = overlaps[i];
        if (!is_contained[it->a_id()] && !is_contained[it->b_id()]) {
            if (overlap_type[i] < 3) {
                continue;
            }

            auto a = nodes_[read_id_to_node_id[it->a_id()] + (it->a_rc() == 0 ? 0 : 1)].get();
            auto _a = a->pair;

            auto b = nodes_[read_id_to_node_id[it->b_id()] + (it->b_rc() == 0 ? 0 : 1)].get();
            auto _b = b->pair;

            if (overlap_type[i] == 3) { // a to b overlap
                Edge* edge = new Edge(edge_id++, overlaps[i], a, b);
                Edge* _edge = new Edge(edge_id++, overlaps[i], _b, _a);

                edge->pair = _edge;
                _edge->pair = edge;

                edges_.push_back(std::unique_ptr<Edge>(edge));
                edges_.push_back(std::unique_ptr<Edge>(_edge));

                a->suffix_edges.push_back(edge);
                _a->prefix_edges.push_back(_edge);
                b->prefix_edges.push_back(edge);
                _b->suffix_edges.push_back(_edge);

            } else if (overlap_type[i] == 4) { // b to a overlap
                Edge* edge = new Edge(edge_id++, overlaps[i], b, a);
                Edge* _edge = new Edge(edge_id++, overlaps[i], _a, _b);

                edge->pair = _edge;
                _edge->pair = edge;

                edges_.push_back(std::unique_ptr<Edge>(edge));
                edges_.push_back(std::unique_ptr<Edge>(_edge));

                b->suffix_edges.push_back(edge);
                _b->prefix_edges.push_back(_edge);
                a->prefix_edges.push_back(edge);
                _a->suffix_edges.push_back(_edge);
            }
        }
    }
}

Graph::~Graph() {
}

void Graph::simplify() {

    // remove nodes wihtout edges
    for (uint32_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i] == nullptr) {
            continue;
        }
        if (nodes_[i]->in_degree() == 0 && nodes_[i]->out_degree() == 0) {
            nodes_[i].reset();
        }
    }

    // remove transitive edges (inspired by Myers 1995 & 2005)
    const double eps = 0.12;
    auto is_similar = [&](double x, double y) -> bool {
        return (x >= y * (1 - eps) && x <= y * (1 + eps));
    };

    std::vector<Edge*> candidate_edge(nodes_.size(), nullptr);

    for (uint32_t v = 0; v < nodes_.size(); v += 2) {
        if (nodes_[v] == nullptr) continue;

        for (const auto& it: nodes_[v]->suffix_edges) {
            candidate_edge[it->end_node->id] = it;
        }

        for (const auto& it: nodes_[v]->suffix_edges) {
            uint32_t w = it->end_node->id;
            for (const auto& it2: nodes_[w]->suffix_edges) {
                uint32_t x = it2->end_node->id;
                if (candidate_edge[x] != nullptr) {
                    fprintf(stderr, "%u (%u) -> %u (%u)-> %u (%u)\n", v, nodes_[v]->length(), w, nodes_[w]->length(), x, nodes_[x]->length());
                    fprintf(stderr, "%u -> %u (%d)\n", v, w, it->length());
                    fprintf(stderr, "%u -> %u (%d)\n", w, x, it2->length());
                    fprintf(stderr, "%u -> %u (%d)\n", v, x, candidate_edge[x]->length());
                    if (is_similar(it->length() + it2->length(), candidate_edge[x]->length())) {
                        fprintf(stderr, "WE DID IT REDIT\n");
                        candidate_edge[x]->mark = true;
                        candidate_edge[x]->pair->mark = true;
                    }
                }
            }

            // remove short overlaps
            /*uint32_t v_length = nodes_[v]->length();
            const double max_ratio = 0.7;
            for (const auto& it2: nodes_[v]->suffix_edges()) {
                if (it == it2) continue;
                if (edges_[it2] == nullptr) continue;
                if (reduce[it] || reduce[it2]) continue;

                if (edges_[it]->length() <= edges_[it2]->length()) {
                    fprintf(stderr, "%u -> %u && %u -> %u\n", v, edges_[it]->end_node_id(), v, edges_[it2]->end_node_id());
                    if ((v_length - edges_[it2]->length()) / (double) (v_length - edges_[it]->length()) < max_ratio) {
                        fprintf(stderr, "WE DID IT REDIT 2nd!\n");
                        reduce[it2] = true;
                    }
                }
            }*/
        }

        for (const auto& it: nodes_[v]->suffix_edges) {
            candidate_edge[it->end_node->id] = nullptr;
        }
    }

    for (const auto& it: nodes_) {
        if (it == nullptr) continue;

        auto it2 = it->prefix_edges.begin();
        while (it2 != it->prefix_edges.end()) {
            if ((*it2)->mark == true) {
                it2 = it->prefix_edges.erase(it2);
            } else {
                ++it2;
            }
        }

        it2 = it->suffix_edges.begin();
        while (it2 != it->suffix_edges.end()) {
            if ((*it2)->mark == true) {
                it2 = it->suffix_edges.erase(it2);
            } else {
                ++it2;
            }
        }
    }

    for (uint32_t i = 0; i < edges_.size(); ++i) {
        if (edges_[i] == nullptr) continue;
        if (edges_[i]->mark == true) {
            edges_[i].reset();
        }
    }
}

void Graph::print() const {

    printf("digraph 1 {\n");
    printf("    overlap = scalexy\n");

    for (const auto& it: nodes_) {
        if (it == nullptr) continue;
        printf("    %d [label = \"%u: +%u,-%u\"", it->id, it->id, it->in_degree(), it->out_degree());
        if (it->rc == true) {
            printf(", style = filled, fillcolor = brown1]\n");
            printf("    %d -> %d [style = dotted, arrowhead = none]\n", it->id, it->id - 1);
        } else {
            printf("]\n");
        }
    }

    for (const auto& it: edges_) {
        if (it == nullptr) continue;
        printf("    %d -> %d [label = \"%d\"]\n", it->begin_node->id, it->end_node->id, it->length());
    }

    printf("}\n");
}

}
