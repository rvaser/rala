/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <set>
#include <list>
#include <stack>
#include <deque>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <set>

#include "read.hpp"
#include "overlap.hpp"
#include "graph.hpp"
#include "preprocess.hpp"

namespace rala {

constexpr double kTransitiveEdgeEps = 0.12;
constexpr uint32_t kMaxBubbleLength = 5000000;
constexpr uint32_t kMinUnitigSize = 6;
constexpr double kMaxOverlapRatio = 0.9;

class Graph::Node {
    public:
        // Node encapsulating read
        Node(uint32_t _id, const std::shared_ptr<Read>& read) :
                id(_id), read_id(read->id()), pair(), sequence(id % 2 == 0 ? read->sequence() : read->rc()),
                prefix_edges(), suffix_edges(), read_ids(1, read->id()), unitig_size(1), mark(false) {

            auto is_unitig = read->name().find("Utg=");
            if (is_unitig != std::string::npos) {
                unitig_size = atoi(read->name().c_str() + is_unitig + 4);
                // fprintf(stderr, "Unitig size = %u\n", unitig_size);
            }
        }
        // Unitig
        Node(uint32_t _id, Node* begin_node, Node* end_node, std::unordered_set<uint32_t>& marked_edges);
        // Circular unitig
        Node(uint32_t _id, Node* begin_node, std::unordered_set<uint32_t>& marked_edges);
        Node(const Node&) = delete;
        const Node& operator=(const Node&) = delete;

        ~Node() {}

        uint32_t length() const {
            return sequence.size();
        }

        uint32_t in_degree() const {
            return prefix_edges.size();
        }

        uint32_t out_degree() const {
            return suffix_edges.size();
        }

        bool is_junction() const {
            return (out_degree() > 1 || in_degree() > 1);
        }

        bool is_tip() const {
            return (out_degree() > 0 && in_degree() == 0 && unitig_size < kMinUnitigSize);
        }

        uint32_t id;
        uint32_t read_id;
        Node* pair;
        std::string sequence;
        std::list<Edge*> prefix_edges;
        std::list<Edge*> suffix_edges;
        std::vector<uint32_t> read_ids;
        uint32_t unitig_size;
        bool mark;
};

class Graph::Edge {
    public:
        Edge(uint32_t _id, const std::shared_ptr<Overlap>& overlap, Node* _begin_node,
            Node* _end_node, uint32_t type) :
                id(_id), pair(), begin_node(_begin_node), end_node(_end_node), length(),
                quality(overlap->quality()), mark(false) {

            uint32_t length_a = id % 2 == 0 ? overlap->a_begin() : overlap->a_length() - overlap->a_end();
            uint32_t length_b = id % 2 == 0 ? overlap->b_begin() : overlap->b_length() - overlap->b_end();

            if (type == 0) { // a to b overlap
                length = length_a - length_b;
            } else { // b to a overlap
                length = length_b - length_a;
            }
        }
        Edge(const Edge&) = delete;
        const Edge& operator=(const Edge&) = delete;

        ~Edge() {}

        std::string label() const {
            return begin_node->sequence.substr(0, length);
        }

        uint32_t matching_bases() const {
            return (quality * (begin_node->length() - length));
        }

        uint32_t id;
        Edge* pair;
        Node* begin_node;
        Node* end_node;
        uint32_t length;
        double quality;
        bool mark;
};

Graph::Node::Node(uint32_t _id, Node* begin_node, Node* end_node, std::unordered_set<uint32_t>& marked_edges) :
        id(_id), read_id(), pair(), sequence(), prefix_edges(), suffix_edges(),
        read_ids(), unitig_size(), mark(false) {

    if (!begin_node->prefix_edges.empty()) {
        begin_node->prefix_edges.front()->end_node = this;
        prefix_edges.push_back(begin_node->prefix_edges.front());
    }

    uint32_t length = 0;
    Node* curr_node = begin_node;
    while (curr_node->id != end_node->id) {
        auto* edge = curr_node->suffix_edges.front();
        edge->mark = true;
        marked_edges.insert(edge->id);

        read_ids.reserve(read_ids.size() + curr_node->read_ids.size());
        read_ids.insert(read_ids.end(), curr_node->read_ids.begin(), curr_node->read_ids.end());

        unitig_size += curr_node->unitig_size;
        length += edge->length;
        sequence += edge->label();

        curr_node->prefix_edges.clear();
        curr_node->suffix_edges.clear();
        curr_node->mark = true;

        curr_node = edge->end_node;
    }

    read_ids.reserve(read_ids.size() + end_node->read_ids.size());
    read_ids.insert(read_ids.end(), end_node->read_ids.begin(), end_node->read_ids.end());

    unitig_size += end_node->unitig_size;
    sequence += end_node->sequence;

    if (!end_node->suffix_edges.empty()) {
        end_node->suffix_edges.front()->begin_node = this;
        end_node->suffix_edges.front()->length += length;
        suffix_edges.push_back(end_node->suffix_edges.front());
    }

    end_node->prefix_edges.clear();
    end_node->suffix_edges.clear();
    end_node->mark = true;
}

Graph::Node::Node(uint32_t _id, Node* begin_node, std::unordered_set<uint32_t>& marked_edges) :
        id(_id), read_id(), pair(), sequence(), prefix_edges(), suffix_edges(),
        read_ids(), unitig_size(), mark(false) {
    // fprintf(stderr, "!!! CIRCULAR UNITIG ALERT !!!\n");

    uint32_t length = 0;
    Node* curr_node = begin_node;
    while (true) {
        // fprintf(stderr, "Curr node = %d\n", curr_node->id);
        auto* edge = curr_node->suffix_edges.front();
        edge->mark = true;
        marked_edges.insert(edge->id);

        read_ids.reserve(read_ids.size() + curr_node->read_ids.size());
        read_ids.insert(read_ids.end(), curr_node->read_ids.begin(), curr_node->read_ids.end());

        unitig_size += curr_node->unitig_size;
        length += edge->length;
        sequence += edge->label();

        curr_node->prefix_edges.clear();
        curr_node->suffix_edges.clear();
        curr_node->mark = true;

        curr_node = edge->end_node;
        if (curr_node->id == begin_node->id) {
            break;
        }
    }
}

std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps) {
    return std::unique_ptr<Graph>(new Graph(reads, overlaps));
}

Graph::Graph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps)
        : nodes_(), edges_(), marked_edges_() {

    fprintf(stderr, "Assembly graph {\n");

    uint64_t max_read_id = 0;
    for (const auto& it: reads) {
        max_read_id = std::max(max_read_id, it->id());
    }

    // create assembly graph
    std::vector<int32_t> read_id_to_node_id(max_read_id + 1, -1);
    uint32_t node_id = 0;
    for (const auto& read: reads) {
        read_id_to_node_id[read->id()] = node_id;

        Node* node = new Node(node_id++, read); // normal read
        Node* _node = new Node(node_id++, read); // reverse complement

        node->pair = _node;
        _node->pair = node;

        nodes_.push_back(std::unique_ptr<Node>(node));
        nodes_.push_back(std::unique_ptr<Node>(_node));
    }

    uint32_t edge_id = 0;
    for (const auto& overlap: overlaps) {

        auto a = nodes_[read_id_to_node_id[overlap->a_id()] + (overlap->a_rc() == 0 ? 0 : 1)].get();
        auto _a = a->pair;

        auto b = nodes_[read_id_to_node_id[overlap->b_id()] + (overlap->b_rc() == 0 ? 0 : 1)].get();
        auto _b = b->pair;

        if (overlap->type() == 3) { // a to b overlap
            Edge* edge = new Edge(edge_id++, overlap, a, b, 0);
            Edge* _edge = new Edge(edge_id++, overlap, _b, _a, 1);

            edge->pair = _edge;
            _edge->pair = edge;

            edges_.push_back(std::unique_ptr<Edge>(edge));
            edges_.push_back(std::unique_ptr<Edge>(_edge));

            a->suffix_edges.push_back(edge);
            _a->prefix_edges.push_back(_edge);
            b->prefix_edges.push_back(edge);
            _b->suffix_edges.push_back(_edge);

        } else if (overlap->type() == 4) { // b to a overlap
            Edge* edge = new Edge(edge_id++, overlap, b, a, 1);
            Edge* _edge = new Edge(edge_id++, overlap, _a, _b, 0);

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

    fprintf(stderr, "  Construction {\n");
    fprintf(stderr, "    number of graph nodes = %zu\n", nodes_.size());
    fprintf(stderr, "    number of graph edges = %zu\n", edges_.size());
    fprintf(stderr, "  }\n");
}

Graph::~Graph() {
    fprintf(stderr, "}\n");
}

void Graph::remove_isolated_nodes() {

    for (auto& node: nodes_) {
        if (node == nullptr) {
            continue;
        }
        if ((node->in_degree() == 0 && node->out_degree() == 0 && node->unitig_size < kMinUnitigSize) || (node->mark == true)) {
            // fprintf(stderr, "Removing isolated node: %d\n", node->id);
            node.reset();
        }
    }
}

void Graph::remove_transitive_edges() {

    fprintf(stderr, "  Transitive edge removal {\n");

    uint32_t num_transitive_edges = 0;
    std::vector<Edge*> candidate_edge(nodes_.size(), nullptr);

    for (const auto& node_x: nodes_) {
        if (node_x == nullptr) continue;

        for (const auto& edge: node_x->suffix_edges) {
            candidate_edge[edge->end_node->id] = edge;
        }

        for (const auto& edge_xy: node_x->suffix_edges) {
            for (const auto& edge_yz: nodes_[edge_xy->end_node->id]->suffix_edges) {
                uint32_t z = edge_yz->end_node->id;
                if (candidate_edge[z] != nullptr && candidate_edge[z]->mark == false) {
                    if (comparable(edge_xy->length + edge_yz->length, candidate_edge[z]->length, kTransitiveEdgeEps)) {
                        candidate_edge[z]->mark = true;
                        candidate_edge[z]->pair->mark = true;
                        marked_edges_.insert(candidate_edge[z]->id);
                        marked_edges_.insert(candidate_edge[z]->pair->id);
                        ++num_transitive_edges;
                    }
                }
            }
        }

        for (const auto& edge: node_x->suffix_edges) {
            candidate_edge[edge->end_node->id] = nullptr;
        }
    }
    remove_marked_edges();

    fprintf(stderr, "    removed %u edges\n", num_transitive_edges);
    fprintf(stderr, "  }\n");
}

void Graph::remove_long_edges() {

    fprintf(stderr, "  Long edge removal {\n");
    uint32_t num_long_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || node->suffix_edges.size() < 2) continue;

        for (const auto& edge1: node->suffix_edges) {
            for (const auto& edge2: node->suffix_edges) {
                if (edge1->id == edge2->id || edge1->mark == true || edge2->mark == true) continue;
                if ((node->length() - edge2->length) / (double) (node->length() - edge1->length) < kMaxOverlapRatio) {
                    edge2->mark = true;
                    edge2->pair->mark = true;
                    marked_edges_.insert(edge2->id);
                    marked_edges_.insert(edge2->pair->id);
                    ++num_long_edges;
                }
            }
        }
    }

    fprintf(stderr, "    removed %u edges\n", num_long_edges);
    fprintf(stderr, "  }\n");

    remove_marked_edges();
}

uint32_t Graph::remove_tips() {

    fprintf(stderr, "  Tip removal {\n");

    uint32_t num_tip_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr) {
            continue;
        }
        // fprintf(stderr, "Considering node %d for tip removal\r", node->id);
        if (!node->is_tip()) {
            continue;
        }

        uint32_t num_removed_edges = 0;

        for (const auto& edge: node->suffix_edges) {
            if (edge->end_node->in_degree() > 1) {
                // fprintf(stderr, "Removing %d\n", edge->begin_node->id);
                edge->mark = true;
                edge->pair->mark = true;
                marked_edges_.insert(edge->id);
                marked_edges_.insert(edge->pair->id);
                ++num_removed_edges;
            }
        }

        if (num_removed_edges == node->suffix_edges.size()) {
            node->mark = true;
            node->pair->mark = true;
        }

        num_tip_edges += num_removed_edges;

        remove_marked_edges();
    }
    remove_isolated_nodes();

    fprintf(stderr, "    removed %u edges\n", num_tip_edges);
    fprintf(stderr, "  }\n");

    return num_tip_edges;
}

void Graph::remove_cycles() {

    fprintf(stderr, "  Cycle removal {\n");

    std::stack<uint32_t> stack;
    std::vector<int32_t> indexes(nodes_.size(), -1);
    std::vector<int32_t> low_links(nodes_.size(), -1);
    std::vector<bool> is_on_stack(nodes_.size(), false);
    int32_t index = 0;

    std::vector<std::vector<uint32_t>> cycles;

    std::function<void(uint32_t)> strong_connect = [&](uint32_t v) -> void {
        indexes[v] = index;
        low_links[v] = index;
        ++index;
        // fprintf(stderr, "Pushing %d\n", v);
        stack.push(v);
        is_on_stack[v] = true;

        for (const auto& edge: nodes_[v]->suffix_edges) {
            uint32_t w = edge->end_node->id;
            if (indexes[w] == -1) {
                strong_connect(w);
                low_links[v] = std::min(low_links[v], low_links[w]);
            } else if (is_on_stack[w]) {
                low_links[v] = std::min(low_links[v], indexes[w]);
            }
        }

        if (low_links[v] == indexes[v]) {
            // new strongly connected component
            std::vector<uint32_t> scc = { v };
            uint32_t w;
            do {
                w = stack.top();
                stack.pop();
                is_on_stack[w] = false;
                scc.push_back(w);
            } while (v != w);

            if (scc.size() > 2) {
                cycles.push_back(scc);
            }
        }
    };

    uint32_t num_cycle_edges = 0;
    do {
        cycles.clear();
        for (const auto& node: nodes_) {
            if (node == nullptr) continue;
            if (indexes[node->id] == -1) {
                strong_connect(node->id);
            }
        }

        // fprintf(stderr, "Number of cycles %zu\n", cycles.size());

        for (const auto& cycle: cycles) {

            Edge* worst_edge = nullptr;
            double min_score = 5;

            for (uint32_t i = 0; i < cycle.size() - 1; ++i) {
                const auto& node = nodes_[cycle[i]];
                for (auto& edge: node->prefix_edges) {
                    if (edge->begin_node->id == cycle[i + 1]) {
                        if (min_score > edge->quality) {
                            min_score = edge->quality;
                            worst_edge = edge;
                        }
                        break;
                    }
                }
            }

            worst_edge->mark = true;
            worst_edge->pair->mark = true;
            marked_edges_.insert(worst_edge->id);
            marked_edges_.insert(worst_edge->pair->id);
            ++num_cycle_edges;
        }

        remove_marked_edges();

    } while (cycles.size() != 0);

    fprintf(stderr, "    removed %u edges\n", num_cycle_edges);
    fprintf(stderr, "  }\n");
}

uint32_t Graph::remove_chimeras() {

    fprintf(stderr, "  Chimera removal {\n");

    std::vector<uint32_t> visited(nodes_.size(), 0);
    uint32_t visited_length = 0;
    std::vector<int32_t> predecessor(nodes_.size(), -1);
    std::deque<uint32_t> queue;

    auto extract_path = [&](std::vector<int32_t>& dst, int32_t sink) -> void {
        int32_t curr_id = sink;
        while (curr_id != -1) {
            dst.push_back(curr_id);
            curr_id = predecessor[curr_id];
        }
    };

    uint32_t num_chimeric_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || node->out_degree() < 2) {
            continue;
        }

        bool found_chim_sink = false;
        int32_t chim_sink = -1, chim_sink_pair = -1;
        int32_t source = node->id;
        std::vector<int32_t> chimeric_edges;

        // DFS
        queue.push_front(source);
        visited[visited_length++] = source;
        while (queue.size() != 0 && !found_chim_sink) {
            int32_t v = queue.front();
            queue.pop_front();
            const auto& curr_node = nodes_[v];

            for (const auto& edge: curr_node->suffix_edges) {
                int32_t w = edge->end_node->id;
                if (predecessor[w] != -1 || w == source) {
                    // Cycle or bubble
                    continue;
                }

                visited[visited_length++] = w;
                predecessor[w] = v;
                queue.push_front(w);

                int32_t w_pair = edge->end_node->pair->id;
                if (predecessor[w_pair] != -1) {
                    // Chimeric link!
                    chim_sink = w;
                    chim_sink_pair = w_pair;
                    found_chim_sink = true;
                    break;
                }
            }
        }

        if (found_chim_sink) {
            // fprintf(stderr, "Source = %d, %d, %d\n", source, chim_sink, chim_sink_pair);

            std::vector<int32_t> path;
            extract_path(path, chim_sink);
            if (path.back() != chim_sink_pair) {
                std::vector<int32_t> other_path;
                extract_path(other_path, chim_sink_pair);

                int32_t ancestor_i = -1, ancestor_j = -1;
                for (uint32_t i = 0; i < path.size(); ++i) {
                    for (uint32_t j = 0; j < other_path.size(); ++j) {
                        if (path[i] == other_path[j]) {
                            ancestor_i = i;
                            ancestor_j = j;
                            break;
                        }
                    }
                    if (ancestor_i != -1) {
                        break;
                    }
                }
                if (ancestor_i != -1) path.resize(ancestor_i + 1);
                std::reverse(path.begin(), path.end());
                if (ancestor_j != -1) other_path.resize(ancestor_j + 1);

                for (uint32_t j = 1; j < other_path.size(); ++j) {
                    path.push_back(nodes_[other_path[j]]->pair->id);
                }
            } else {
                std::reverse(path.begin(), path.end());
            }

            std::vector<uint32_t> chimeric_edges;
            find_removable_edges(chimeric_edges, path, true);

            if (chimeric_edges.size() > 0) {
                // remove chimeric edges
                // for (const auto& pid: path) fprintf(stderr, "%d -> ", pid);
                // fprintf(stderr, "\n");
                for (const auto& edge_id: chimeric_edges) {
                    // fprintf(stderr, "%d\n", edge_id);
                    // fprintf(stderr, "Removing: %d -> %d\n", edges_[edge_id]->begin_node->id, edges_[edge_id]->end_node->id);
                    edges_[edge_id]->mark = true;
                    edges_[edge_id]->pair->mark = true;
                    marked_edges_.insert(edge_id);
                    marked_edges_.insert(edges_[edge_id]->pair->id);
                }
                ++num_chimeric_edges;
                remove_marked_edges();
            } else {
                // for (const auto& pid: path) fprintf(stderr, "%d -> ", pid);
                // fprintf(stderr, "\n");
            }
        }

        queue.clear();
        for (uint32_t i = 0; i < visited_length; ++i) {
            predecessor[visited[i]] = -1;
        }
        visited_length = 0;
    }

    fprintf(stderr, "    removed %u edges\n", num_chimeric_edges);
    fprintf(stderr, "  }\n");

    return num_chimeric_edges;
}

uint32_t Graph::remove_bubbles() {

    fprintf(stderr, "  Bubble removal {\n");

    std::vector<uint32_t> distance(nodes_.size(), 0);
    std::vector<uint32_t> visited(nodes_.size(), 0);
    uint32_t visited_length = 0;
    std::vector<int32_t> predecessor(nodes_.size(), -1);
    std::deque<uint32_t> queue;

    auto extract_path = [&](std::vector<int32_t>& dst, int32_t source, int32_t sink) -> void {
        int32_t curr_id = sink;
        while (curr_id != source) {
            dst.push_back(curr_id);
            curr_id = predecessor[curr_id];
        }
        dst.push_back(source);
        std::reverse(dst.begin(), dst.end());
    };

    auto calculate_path_length = [&](const std::vector<int32_t>& path) -> uint32_t {
        uint32_t path_length = nodes_[path.back()]->length();
        for (uint32_t i = 0; i < path.size() - 1; ++i) {
            for (const auto& edge: nodes_[path[i]]->suffix_edges) {
                if (edge->end_node->id == (uint32_t) path[i + 1]) {
                    path_length += edge->length;
                    break;
                }
            }
        }
        return path_length;
    };

    auto is_valid_bubble = [&](const std::vector<int32_t>& path, const std::vector<int32_t>& other_path) -> bool {
        std::set<int32_t> node_set;
        for (const auto& id: path) node_set.insert(id);
        for (const auto& id: other_path) node_set.insert(id);
        if (path.size() + other_path.size() - 2 != node_set.size()) {
            return false;
        }
        for (const auto& id: path) {
            uint32_t pair_id = (id % 2 == 0) ? id + 1 : id - 1;
            if (node_set.count(pair_id) != 0) {
                return false;
            }
        }
        uint32_t path_length = calculate_path_length(path);
        uint32_t other_path_length = calculate_path_length(other_path);
        if (std::min(path_length, other_path_length) / (double) std::max(path_length, other_path_length) < 0.8) {
            for (uint32_t i = 1; i < other_path.size() - 1; ++i) {
                if (nodes_[other_path[i]]->in_degree() > 1 || nodes_[other_path[i]]->out_degree() > 1) {
                    return false;
                }
            }
            for (uint32_t i = 1; i < path.size() - 1; ++i) {
                if (nodes_[path[i]]->in_degree() > 1 || nodes_[path[i]]->out_degree() > 1) {
                    return false;
                }
            }
        }
        return true;
    };

    uint32_t num_bubbles_popped = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr || node->out_degree() < 2) continue;

        bool found_sink = false;
        int32_t sink = 0, sink_other_predecesor = 0;
        int32_t source = node->id;

        // BFS
        queue.push_back(source);
        visited[visited_length++] = source;
        while (queue.size() != 0 && !found_sink) {
            int32_t v = queue.front();
            const auto& curr_node = nodes_[v];

            queue.pop_front();

            for (const auto& edge: curr_node->suffix_edges) {
                int32_t w = edge->end_node->id;

                if (w == source) {
                    // Cycle
                    continue;
                }

                if (distance[v] + edge->length > kMaxBubbleLength) {
                    // Out of reach
                    continue;
                }

                distance[w] = distance[v] + edge->length;
                visited[visited_length++] = w;
                queue.push_back(w);

                if (predecessor[w] != -1) {
                    sink = w;
                    sink_other_predecesor = v;
                    found_sink = true;
                    break;
                }

                predecessor[w] = v;
            }
        }

        if (found_sink) {
            // fprintf(stderr, "Source = %u, sink = %u, sink_predecesors = [%u, %u]\n", source, sink, predecessor[sink], sink_other_predecesor);

            std::vector<int32_t> path, other_path;
            extract_path(path, source, sink);
            other_path.push_back(sink);
            extract_path(other_path, source, sink_other_predecesor);

            /*fprintf(stderr, "Path 1:");
            for (const auto& it: path) fprintf(stderr, " %d", it);
            fprintf(stderr, "\n");

            fprintf(stderr, "Path 2:");
            for (const auto& it: other_path) fprintf(stderr, " %d", it);
            fprintf(stderr, "\n");*/

            if (!is_valid_bubble(path, other_path)) {
                // fprintf(stderr, "Not valid bubble!\n");
            } else {
                uint32_t path_num_reads = 0;
                for (const auto& it: path) path_num_reads += nodes_[it]->unitig_size;

                uint32_t other_path_num_reads = 0;
                for (const auto& it: other_path) other_path_num_reads += nodes_[it]->unitig_size;

                std::vector<uint32_t> edges_for_removal;
                if (path_num_reads > other_path_num_reads) {
                    // fprintf(stderr, "2\n");
                    find_removable_edges(edges_for_removal, other_path);
                } else {
                    // fprintf(stderr, "1\n");
                    find_removable_edges(edges_for_removal, path);
                }

                for (const auto& edge_id: edges_for_removal) {
                    // fprintf(stderr, "%d\n", edge_id);
                    // fprintf(stderr, "Removing: %d -> %d\n", edges_[edge_id]->begin_node->id, edges_[edge_id]->end_node->id);
                    edges_[edge_id]->mark = true;
                    edges_[edge_id]->pair->mark = true;
                    marked_edges_.insert(edge_id);
                    marked_edges_.insert(edges_[edge_id]->pair->id);
                }
                if (!edges_for_removal.empty()) {
                    remove_marked_edges();
                    ++num_bubbles_popped;
                }
            }
        }

        queue.clear();
        for (uint32_t i = 0; i < visited_length; ++i) {
            distance[visited[i]] = 0;
            predecessor[visited[i]] = -1;
        }
        visited_length = 0;
    }

    remove_isolated_nodes();

    fprintf(stderr, "    popped %d bubbles\n", num_bubbles_popped);
    fprintf(stderr, "  }\n");

    return num_bubbles_popped;
}

uint32_t Graph::create_unitigs() {

    fprintf(stderr, "  Creating unitigs {\n");

    uint32_t node_id = nodes_.size();
    std::vector<bool> visited(nodes_.size(), false);
    std::vector<std::unique_ptr<Node>> new_nodes;

    uint32_t num_unitigs_created = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || visited[node->id] || node->is_junction()) continue;
        // fprintf(stderr, "Considering node %d for unittiging\r", node->id);

        bool is_circular = false;
        auto bnode = node.get();
        while (!bnode->is_junction()) {
            visited[bnode->id] = true;
            visited[bnode->pair->id] = true;
            if (bnode->in_degree() == 0 || bnode->prefix_edges.front()->begin_node->is_junction()) {
                break;
            }
            bnode = bnode->prefix_edges.front()->begin_node;
            if (bnode->id == node->id) {
                is_circular = true;
                break;
            }
        }

        auto enode = node.get();
        while (!enode->is_junction()) {
            visited[enode->id] = true;
            visited[enode->pair->id] = true;
            if (enode->out_degree() == 0 || enode->suffix_edges.front()->end_node->is_junction()) {
                break;
            }
            //++unitig_size;
            enode = enode->suffix_edges.front()->end_node;
            if (enode->id == node->id) {
                is_circular = true;
                break;
            }
        }

        // normal
        Node* unitig = nullptr;
        // reverse_complement
        Node* _unitig = nullptr;

        if (is_circular) {
            unitig = new Node(node_id++, bnode, marked_edges_);
            _unitig = new Node(node_id++, bnode->pair, marked_edges_);
        } else if (bnode->id != enode->id) {
            unitig = new Node(node_id++, bnode, enode, marked_edges_);
            _unitig = new Node(node_id++, enode->pair, bnode->pair, marked_edges_);
        }

        if (unitig != nullptr && _unitig != nullptr) {
            unitig->pair = _unitig;
            _unitig->pair = unitig;

            new_nodes.push_back(std::unique_ptr<Node>(unitig));
            new_nodes.push_back(std::unique_ptr<Node>(_unitig));

            // fprintf(stderr, "Unitig: %d -> %d && %d -> %d\n", bnode->id, enode->id, enode->pair->id, bnode->pair->id);
            ++num_unitigs_created;
        }
    }

    for (auto& node: new_nodes) {
        nodes_.push_back(std::move(node));
    }

    remove_marked_edges();
    remove_isolated_nodes();

    fprintf(stderr, "    created %u new unitigs\n", num_unitigs_created);
    fprintf(stderr, "  }\n");

    return num_unitigs_created;
}

void Graph::print_contigs() const {

    fprintf(stderr, "  Contig creation {\n");

    uint32_t contig_id = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr || node->id % 2 == 0 || node->unitig_size < kMinUnitigSize) continue;
        fprintf(stderr, "    >Contig_%d_(Utg:%u), length = %zu (%d -> %d)\n", contig_id,
            node->unitig_size, node->sequence.size(), node->read_ids.front(),
            node->read_ids.back());
        fprintf(stdout, ">Contig_%u_(Utg=%u:Len=%lu)\n%s\n", contig_id++, node->unitig_size, node->sequence.size(), node->sequence.c_str());
    }

    fprintf(stderr, "  }\n");
}

int32_t Graph::find_edge(uint32_t src, uint32_t dst) {
    for (const auto& edge: nodes_[src]->suffix_edges) {
        if (edge->end_node->id == dst) {
            return edge->id;
        }
    }
    return -1;
}

void Graph::find_removable_edges(std::vector<uint32_t>& dst, const std::vector<int32_t>& path,
    bool chimeric) {

    // find first node with multiple in edges
    int32_t pref = -1;
    for (uint32_t i = 1; i < path.size() - 1; ++i) {
        if (nodes_[path[i]]->in_degree() > 1) {
            pref = i;
            break;
        }
    }
    // find last node with multiple out edges
    int32_t suff = -1;
    for (uint32_t i = 1; i < path.size() - 1; ++i) {
        if (nodes_[path[i]]->out_degree() > 1) {
            suff = i;
        }
    }

    if (pref == -1 && suff == -1) {
        if (chimeric) {
            // remove first or last edge
            const auto& first_edge = edges_[find_edge(path[0], path[1])];
            const auto& last_edge = edges_[find_edge(path[path.size() - 2], path[path.size() - 1])];
            fprintf(stderr, "%d %d\n", first_edge->id, first_edge->begin_node->length() - first_edge->length);
            fprintf(stderr, "%d %d\n", last_edge->id, last_edge->begin_node->length() - last_edge->length);
            if ((first_edge->begin_node->length() - first_edge->length) > (last_edge->begin_node->length() - last_edge->length)) {
                dst.push_back(last_edge->id);
            } else {
                dst.push_back(first_edge->id);
            }
        } else {
            // remove whole path
            for (uint32_t i = 0; i < path.size() - 1; ++i) {
                dst.push_back(find_edge(path[i], path[i + 1]));
            }
        }
        return;
    }

    if (pref != -1 && nodes_[path[pref]]->out_degree() > 1) return;
    if (suff != -1 && nodes_[path[suff]]->in_degree() > 1) return;

    if (pref == -1) {
        // remove everything after last suff node
        for (uint32_t i = suff; i < path.size() - 1; ++i) {
            dst.push_back(find_edge(path[i], path[i + 1]));
        }
    } else if (suff == -1) {
        // remove everything before first pref node
        for (int32_t i = 0; i < pref; ++i) {
            dst.push_back(find_edge(path[i], path[i + 1]));
        }
    } else if (suff < pref) {
        // remove everything between last suff and first pref node
        for (int32_t i = suff; i < pref; ++i) {
            dst.push_back(find_edge(path[i], path[i + 1]));
        }
    } else if (chimeric && suff >= pref && nodes_[path[0]]->in_degree() == 0) {
        // remove everything after last suff node
        for (uint32_t i = suff; i < path.size() - 1; ++i) {
            dst.push_back(find_edge(path[i], path[i + 1]));
        }
        // remove everything before first pref node
        for (int32_t i = 0; i < pref; ++i) {
            dst.push_back(find_edge(path[i], path[i + 1]));
        }
    }
}

void Graph::remove_marked_edges() {

    auto delete_edges = [&](std::list<Edge*>& edges) -> void {
        auto edge = edges.begin();
        while (edge != edges.end()) {
            if ((*edge)->mark == true) {
                edge = edges.erase(edge);
            } else {
                ++edge;
            }
        }
    };

    std::unordered_set<uint32_t> marked_nodes;
    for (const auto& it: marked_edges_) {
        marked_nodes.insert(edges_[it]->begin_node->id);
        marked_nodes.insert(edges_[it]->end_node->id);
    }

    for (const auto& it: marked_nodes) {
        delete_edges(nodes_[it]->prefix_edges);
        delete_edges(nodes_[it]->suffix_edges);
    }

    for (const auto& it: marked_edges_) {
        edges_[it].reset();
    }
    marked_edges_.clear();

    /*for (const auto& node: nodes_) {
        if (node == nullptr) continue;
        delete_edges(node->prefix_edges);
        delete_edges(node->suffix_edges);
    }
    for (auto& edge: edges_) {
        if (edge == nullptr) continue;
        if (edge->mark == true) {
            edge.reset();
        }
    }*/
}

void Graph::print_csv(std::string path, const std::vector<std::shared_ptr<ReadInfo>>& read_infos) const {

    auto graph_file = fopen(path.c_str(), "w");

    for (const auto& node: nodes_) {
        if (node == nullptr || node->id % 2 == 0) continue;
        fprintf(graph_file, "%u [%u] {%d} U:%d,%u [%u] {%d} U:%d,0,-\n",
            node->id, node->length(), node->read_id, node->unitig_size,
            node->pair->id, node->pair->length(), node->pair->read_id, node->pair->unitig_size);
    }

    for (const auto& edge: edges_) {
        if (edge == nullptr) continue;
        fprintf(graph_file, "%u [%u] {%d} U:%d,%u [%u] {%d} U:%d,1,%d %d %g\n",
            edge->begin_node->id, edge->begin_node->length(), edge->begin_node->read_id, edge->begin_node->unitig_size,
            edge->end_node->id, edge->end_node->length(), edge->end_node->read_id, edge->end_node->unitig_size,
            edge->id, edge->length, edge->quality);
    }

    fclose(graph_file);
}

void Graph::print_knots(const std::vector<std::shared_ptr<ReadInfo>>& read_infos, double median) const {

    std::vector<bool> visited(edges_.size(), false);

    for (const auto& node: nodes_) {
        if (node == nullptr || node->suffix_edges.size() < 2) continue;
        if (read_infos[node->read_id] == nullptr) continue;

        std::vector<uint16_t> graph1(read_infos[node->read_id]->coverage_graph());
        uint32_t begin1 = read_infos[node->read_id]->begin();
        uint32_t end1 = read_infos[node->read_id]->end();
        if (node->id % 2 != 0) {
            std::reverse(graph1.begin(), graph1.end());
            uint32_t tmp = begin1;
            begin1 = graph1.size() - end1;
            end1 = graph1.size() - tmp;
        }

        for (const auto& edge: node->suffix_edges) {
            if (visited[edge->id]) {
                continue;
            }
            visited[edge->id] = true;
            visited[edge->pair->id] = true;

            uint32_t id = edge->end_node->read_id;
            if (read_infos[id] == nullptr) {
                continue;
            }

            std::vector<uint16_t> graph2(read_infos[id]->coverage_graph());
            uint32_t begin2 = read_infos[id]->begin();
            uint32_t end2 = read_infos[id]->end();
            if (edge->end_node->id % 2 != 0) {
                std::reverse(graph2.begin(), graph2.end());
                uint32_t tmp = begin2;
                begin2 = graph2.size() - end2;
                end2 = graph2.size() - tmp;
            }

            std::ofstream out("graphs/e" + std::to_string(edge->id));
            out << "x " << node->read_id << " " << id << " median diff" << std::endl;
            for (uint32_t i = 0; i < begin1 + edge->length + (end2 - begin2); ++i) {
                uint32_t g1 = (i < graph1.size() ? graph1[i] : 0 ), g2 = (i < begin1 + edge->length ? 0 : graph2[i - (edge->length + begin1) + begin2]);
                out << i << " " << g1 << " " << g2 << " " << median << " " << (int32_t) (g1 - g2) << std::endl;
            }
            out.close();
        }
    }
}

void Graph::remove_selected_nodes_and_edges() {

    uint32_t num_chimeras = 0;
    std::set<uint32_t> selected_nodes = {};
    for (const auto& node: nodes_) {
        if (node == nullptr) continue;
        if (selected_nodes.count(node->id) != 0) {
            node->mark = true;
            node->pair->mark = true;
            for (const auto& edge: node->suffix_edges) {
                edge->mark = true;
                edge->pair->mark = true;
                marked_edges_.insert(edge->id);
                marked_edges_.insert(edge->pair->id);
            }
            for (const auto& edge: node->prefix_edges) {
                edge->mark = true;
                edge->pair->mark = true;
                marked_edges_.insert(edge->id);
                marked_edges_.insert(edge->pair->id);
            }

            ++num_chimeras;
        }
    }

    fprintf(stderr, "Num selected nodes = %u\n", num_chimeras);

    std::set<uint32_t> selected_edges = {};
    for (const auto& edge: edges_) {
        if (edge == nullptr) continue;
        if (selected_edges.count(edge->id) != 0) {
            edge->mark = true;
            edge->pair->mark = true;
            marked_edges_.insert(edge->id);
            marked_edges_.insert(edge->pair->id);
        }
    }

    remove_marked_edges();
    remove_isolated_nodes();
}

}
