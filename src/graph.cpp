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

#include "read.hpp"
#include "overlap.hpp"
#include "graph.hpp"

namespace RALA {

constexpr uint32_t kMinOverlapLength = 2000;
constexpr uint32_t kMinMatchingBases = 100;
constexpr double kMinMatchingBasesPerc = 0.055;
constexpr double kChimericRatio = 1.85;
constexpr double kMinMatchingBasesRatio = 2.5;
constexpr double kOverlapQualityRatio = 2.57;
constexpr double kOverlapLengthRatio = 5.17;
constexpr uint32_t kMinCoverage = 2; // 3;
constexpr uint32_t kMaxOverhang = 1000;
constexpr double kMaxOverhangToOverlapRatio = 0.8;
constexpr double kTransitiveEdgeEps = 0.12;
constexpr uint32_t kMaxBubbleLength = 50000;
constexpr uint32_t kMinUnitigSize = 5;

bool isSimilar(double a, double b, double eps) {
    return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
        (b >= a * (1 - eps) && b <= a * (1 + eps));
};

template<class T>
void shrinkVector(std::vector<std::shared_ptr<T>>& dst, uint64_t start, const std::vector<bool>& is_valid) {

    uint64_t i = start, j = start;
    for (; i < dst.size(); ++i) {
        if (is_valid[dst[i]->id()]) {
            continue;
        }

        j = std::max(j, i);
        while (j < dst.size() && !is_valid[dst[j]->id()]) {
            ++j;
        }

        if (j >= dst.size()) {
            break;
        } else if (i != j) {
            dst[i].swap(dst[j]);
        }
    }
    if (i < dst.size()) {
        dst.resize(i);
    }
}

uint32_t classifyOverlap(const std::shared_ptr<Overlap>& overlap) {

    uint32_t left_overhang = std::min(overlap->a_begin(), overlap->b_begin());
    uint32_t right_overhang = std::min(overlap->a_length() - overlap->a_end(),
        overlap->b_length() - overlap->b_end());

    uint32_t a_len = overlap->a_end() - overlap->a_begin();
    uint32_t b_len = overlap->b_end() - overlap->b_begin();

    if (left_overhang > kMaxOverhang * 1.5 || right_overhang > kMaxOverhang * 1.5 ||
        a_len < (a_len + left_overhang + right_overhang) * kMaxOverhangToOverlapRatio ||
        b_len < (b_len + left_overhang + right_overhang) * kMaxOverhangToOverlapRatio) {
        return 0; // internal match
    }
    if (overlap->a_begin() <= overlap->b_begin() && (overlap->a_length() - overlap->a_end()) <= (overlap->b_length() - overlap->b_end())) {
        return 1; // a contained
    }
    if (overlap->a_begin() >= overlap->b_begin() && (overlap->a_length() - overlap->a_end()) >= (overlap->b_length() - overlap->b_end())) {
        return 2; // b contained
    }
    if (overlap->a_begin() > overlap->b_begin()) {
        return 3; // a to b overlap
    }

    return 4; // b to a overlap
}

// call after hit sorting!!
void printCoverageGraph(const std::vector<uint32_t>& hits, const std::string& path) {

    std::vector<uint32_t> coverage_graph((hits.back() >> 1) + 1, 0);
    int32_t coverage = 0;
    uint32_t last = 0;
    for (const auto& hit: hits) {
        if (coverage > 0) {
            for (uint32_t i = last; i < (hit >> 1); ++i) {
                coverage_graph[i] += coverage;
            }
            last = hit >> 1;
        }
        if (hit & 1) --coverage;
        else ++coverage;
    }

    std::ofstream out(path);
    //out.open(path);
    for (uint32_t i = 0; i < coverage_graph.size(); ++i) {
        out << i << " " << coverage_graph[i] << std::endl;
    }
    out.close();
}

// call before hit sorting!!
void printPileGraph(const std::vector<uint32_t>& hits, const std::string& path) {

    std::ofstream out(path);
    // out.open(path);
    for (uint32_t i = 0; i < hits.size(); i+= 2) {
        out << (hits[i] >> 1) << " " << (hits[i+1] >> 1) << std::endl;
    }
    out.close();
};

void prefilterData(std::vector<bool>& is_valid_read, std::vector<bool>& is_valid_overlap,
    const std::string& overlaps_path, uint32_t overlap_type) {

    uint32_t size = 1024 * 1024 * 1024; // ~ 1GB
    auto reader = overlap_type == 0 ?
        BIOPARSER::createReader<Overlap, BIOPARSER::MhapReader>(overlaps_path) :
        BIOPARSER::createReader<Overlap, BIOPARSER::PafReader>(overlaps_path);
    std::vector<std::unique_ptr<Overlap>> overlaps;

    while (true) {
        auto status = reader->read_objects(overlaps, size);

        is_valid_overlap.resize(is_valid_overlap.size() + overlaps.size(), true);

        for (const auto& it: overlaps) {
            uint32_t max_id = std::max(it->a_id(), it->b_id());
            if (is_valid_read.size() <= max_id) {
                is_valid_read.resize(max_id + 1, true);
            }

            if (it->a_end() - it->a_begin() < kMinOverlapLength ||
                it->b_end() - it->b_begin() < kMinOverlapLength ||
                it->matching_bases() < kMinMatchingBases) {
                is_valid_overlap[it->id()] = false;
                continue;
            }

            if (it->a_length() >> 1 > it->b_length()) {
                if (it->b_begin() > kMaxOverhang >> 2 || (it->b_length() - it->b_end()) > kMaxOverhang >> 2 || it->b_end() - it->b_begin() < it->b_length() * kMaxOverhangToOverlapRatio) {
                    continue;
                }
                if (it->a_begin() - it->b_begin() > kMaxOverhang << 1 && (it->a_length() - it->a_end() - (it->b_length() - it->b_end())) > kMaxOverhang << 1) {
                    is_valid_read[it->b_id()] = false;
                }
            } else if (it->a_length() < it->b_length() >> 1 ){
                if (it->a_begin() > kMaxOverhang >> 2 || (it->a_length() - it->a_end()) > kMaxOverhang >> 2 || it->a_end() - it->a_begin() < it->a_length() * kMaxOverhangToOverlapRatio) {
                    continue;
                }
                if (it->b_begin() - it->a_begin() > kMaxOverhang << 1 && (it->b_length() - it->b_end() - (it->a_length() - it->a_end())) > kMaxOverhang << 1) {
                    is_valid_read[it->a_id()] = false;
                }
            }

            /*if (it->a_begin() <= it->b_begin() && (it->a_length() - it->a_end()) <= (it->b_length() - it->b_end())) {
                is_valid_read[it->a_id()] = false;
            }
            else if (it->a_begin() <= it->b_begin() && (it->a_length() - it->a_end()) <= (it->b_length() - it->b_end())) {
                is_valid_read[it->b_id()] = false;
            }*/
        }

        overlaps.clear();

        if (status == false) {
            break;
        }
    }

    uint64_t rtot = 0, otot = 0;
    for (const auto& it: is_valid_read) if (it == true) ++rtot;
    for (const auto& it: is_valid_overlap) if (it == true) ++otot;

    fprintf(stderr, "Valid reads = %lu / %lu\n", rtot, is_valid_read.size());
    fprintf(stderr, "valid overlaps = %lu / %lu\n", otot, is_valid_overlap.size());
}

void preprocessData(std::vector<std::shared_ptr<Read>>& reads, std::vector<std::shared_ptr<Overlap>>& overlaps,
    std::vector<bool>& is_valid_read, std::vector<bool>& is_valid_overlap, const std::string& reads_path,
    const std::string& overlaps_path, uint32_t overlap_type) {

    uint32_t size = 1024 * 1024 * 1024; // ~ 1GB
    auto reader = overlap_type == 0 ?
        BIOPARSER::createReader<Overlap, BIOPARSER::MhapReader>(overlaps_path) :
        BIOPARSER::createReader<Overlap, BIOPARSER::PafReader>(overlaps_path);

    std::vector<std::vector<uint32_t>> hits(is_valid_read.size());

    while (true) {
        auto status = reader->read_objects(overlaps, size);

        for (const auto& it: overlaps) {
            if (!is_valid_overlap[it->id()]) continue;
            if (!is_valid_read[it->a_id()] || !is_valid_read[it->b_id()]) {
                is_valid_overlap[it->id()] = false;
                continue;
            }

            uint32_t begin = it->a_rc() ? it->a_length() - it->a_end() : it->a_begin();
            uint32_t end = it->a_rc() ? it->a_length() - it->a_begin() : it->a_end();
            hits[it->a_id()].push_back(begin << 1 | 0);
            hits[it->a_id()].push_back(end << 1 | 1);

            begin = it->b_rc() ? it->b_length() - it->b_end() : it->b_begin();
            end = it->b_rc() ? it->b_length() - it->b_begin() : it->b_end();
            hits[it->b_id()].push_back(begin << 1 | 0);
            hits[it->b_id()].push_back(end << 1 | 1);
        }

        overlaps.clear();

        if (status == false) {
            break;
        }
    }
    reader.reset();

    std::vector<std::pair<uint32_t, uint32_t>> regions(is_valid_read.size());
    for (uint32_t i = 0; i < hits.size(); ++i) {
        if (hits[i].empty()) {
            is_valid_read[i] = false;
            continue;
        }
        std::sort(hits[i].begin(), hits[i].end());

        int32_t coverage = 0, min_coverage = kMinCoverage;
        int32_t begin = 0, max_begin = 0, max_end = 0;
        for (const auto& hit: hits[i]) {
            int32_t old_coverage = coverage;
            if (hit & 1) {
                --coverage;
            } else {
                ++coverage;
            }
            if (old_coverage < min_coverage && coverage >= min_coverage) {
                begin = hit >> 1;
            } else if (old_coverage >= min_coverage && coverage < min_coverage) {
                int32_t end = (hit >> 1);
                int32_t length = end - begin;
                if (length > max_end - max_begin) {
                    max_begin = begin;
                    max_end = end;
                }
            }
        }

        if (max_end - max_begin > 0) {
            regions[i] = std::make_pair(max_begin, max_end);
        } else {
            is_valid_read[i] = false;
        }
    }

    // 2nd run - chimeric test
    uint32_t chim = 0, brep = 0;
    for (uint32_t i = 0; i < hits.size(); ++i) {
        if (is_valid_read[i] == false) continue;

        int32_t begin = regions[i].first;
        int32_t length = regions[i].second - regions[i].first;
        //int32_t begin = 0;
        //int32_t length = reads[i]->sequence().size();

        int32_t left_border = begin + length / 4;
        int32_t right_border = begin + 3 * length / 4;

        int32_t left_check = begin + length * 0.05;
        bool lcheck = true;
        int32_t left_check_sum = 0, left_check_last = 0;
        int32_t right_check = begin + length * 0.95;
        int32_t right_check_sum = 0, right_check_last = right_check;

        int32_t left_max = 0;
        int32_t right_max = 0;
        int32_t middle_min = 10000;

        int32_t coverage = 0;
        for (const auto& hit: hits[i]) {
            if (hit & 1) {
                --coverage;
            } else {
                ++coverage;
            }

            int32_t pos = hit >> 1;

            if (pos < left_check) {
                left_check_sum += coverage * (pos - left_check_last);
                left_check_last = pos;
            }
            if (pos > left_check && lcheck) {
                lcheck = false;
                left_check_sum += coverage * (left_check - left_check_last);
            }

            if (pos < left_border) {
                left_max = std::max(left_max, coverage);
            }
            if (pos > left_border && pos < right_border) {
                middle_min = std::min(middle_min, coverage);
            }
            if (pos > right_border) {
                right_max = std::max(right_max, coverage);
            }

            if (pos > right_check) {
                right_check_sum += coverage * (pos - right_check_last);
                right_check_last = pos;
            }
        }
        right_check_sum += coverage * (right_check - right_check_last);

        if (left_max / (double) middle_min > kChimericRatio && right_max / (double) middle_min > kChimericRatio) {
            chim++;
            is_valid_read[i] = false;
            // printCoverageGraph(hits[i], "graphs/c" + std::to_string(i));
        } else if (left_max / (double) right_max > kChimericRatio && left_max / (left_check_sum / (double) (left_check - begin)) < kChimericRatio) {
            brep++;
            // is_valid_read[i] = false:
            //printCoverageGraph(hits[i], "graphs/l" + std::to_string(i));
        } else if (right_max / (double) left_max > kChimericRatio && right_max / (right_check_sum / (double) (left_check - begin)) < kChimericRatio) {
            brep++;
            // is_valid_read[i] = false:
            //printCoverageGraph(hits[i], "graphs/r" + std::to_string(i));
        }

        std::vector<uint32_t>().swap(hits[i]);
    }
    std::vector<std::vector<uint32_t>>().swap(hits);

    fprintf(stderr, "Removed %d chimeric reads\n", chim);
    fprintf(stderr, "Removed %d unbridged read repeats\n", brep);

    // reading valid overlaps into memory
    auto oreader = overlap_type == 0 ?
        BIOPARSER::createReader<Overlap, BIOPARSER::MhapReader>(overlaps_path) :
        BIOPARSER::createReader<Overlap, BIOPARSER::PafReader>(overlaps_path);

    uint64_t votot = 0;
    while (true) {
        uint64_t start = overlaps.size();
        auto status = oreader->read_objects(overlaps, size);

        for (auto& it: overlaps) {
            if (!is_valid_overlap[it->id()]) continue;
            if (!is_valid_read[it->a_id()] || !is_valid_read[it->b_id()]) {
                is_valid_overlap[it->id()] = false;
                continue;
            }

            bool is_valid = it->update(
                regions[it->a_id()].first,
                regions[it->a_id()].second,
                regions[it->b_id()].first,
                regions[it->b_id()].second
            );

            if (!is_valid ||
                it->a_end() - it->a_begin() < kMinOverlapLength ||
                it->b_end() - it->b_begin() < kMinOverlapLength ||
                it->matching_bases() < kMinMatchingBases ||
                it->quality() < kMinMatchingBasesPerc) {

                is_valid_overlap[it->id()] = false;
                continue;
            }

            it->set_type(classifyOverlap(it));
            switch (it->type()) {
                case 0:
                    is_valid_overlap[it->id()] = false;
                    break;
                case 1:
                    is_valid_read[it->a_id()] = false;
                    is_valid_overlap[it->id()] = false;
                    break;
                case 2:
                    is_valid_read[it->b_id()] = false;
                    is_valid_overlap[it->id()] = false;
                    break;
                default:
                    ++votot;
                    break;
            }
        }

        shrinkVector(overlaps, start, is_valid_overlap);

        if (status == false) {
            break;
        }
    }
    oreader.reset();

    fprintf(stderr, "Valid overlaps after preprocessing = %lu\n", votot);

    // check if all non valid overlaps are deleted
    bool shrink = false;
    for (const auto& it: overlaps) {
        if (!is_valid_read[it->a_id()] || !is_valid_read[it->b_id()]) {
            shrink = true;
            is_valid_overlap[it->id()] = false;
            --votot;
        }
    }

    if (shrink) {
        shrinkVector(overlaps, 0, is_valid_overlap);
    }

    fprintf(stderr, "Valid overlaps after preprocessing = %lu\n", votot);

    auto rreader = BIOPARSER::createReader<Read, BIOPARSER::FastqReader>(reads_path);
    uint64_t vrtot = 0;

    while (true) {
        uint64_t start = reads.size();
        auto status = rreader->read_objects(reads, size);

        for (auto& it: reads) {
            if (!is_valid_read[it->id()]) {
                continue;
            }
            it->trim_sequence(regions[it->id()].first, regions[it->id()].second);
            ++vrtot;
        }

        shrinkVector(reads, start, is_valid_read);

        if (status == false) {
            break;
        }
    }
    rreader.reset();
    fprintf(stderr, "Reads size = %zu, overlaps size = %zu\n", reads.size(), overlaps.size());

    fprintf(stderr, "Valid reads after preprocessing = %lu\n", vrtot);
}

class Graph::Node {
    public:
        // Node encapsulating read
        Node(uint32_t _id, const std::shared_ptr<Read>& read) :
                id(_id), read_id(read->id()), pair(), sequence(id % 2 == 0 ? read->sequence() : read->rc()),
                prefix_edges(), suffix_edges(), unitig_size(1), mark(false) {
        }
        // Unitig
        Node(uint32_t _id, Node* begin_node, Node* end_node);
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

Graph::Node::Node(uint32_t _id, Node* begin_node, Node* end_node) :
        id(_id), read_id(), pair(), sequence(), prefix_edges(), suffix_edges(),
        unitig_size(), mark(false) {

    if (!begin_node->prefix_edges.empty()) {
        begin_node->prefix_edges.front()->end_node = this;
        prefix_edges.push_back(begin_node->prefix_edges.front());
    }

    uint32_t length = 0;
    Node* curr_node = begin_node;
    while (curr_node->id != end_node->id) {
        auto* edge = curr_node->suffix_edges.front();
        edge->mark = true;

        unitig_size += curr_node->unitig_size;
        length += edge->length;
        sequence += edge->label();

        curr_node->prefix_edges.clear();
        curr_node->suffix_edges.clear();
        curr_node->mark = true;

        curr_node = edge->end_node;
    }

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

std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps) {
    return std::unique_ptr<Graph>(new Graph(reads, overlaps));
}

Graph::Graph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps)
        : nodes_(), edges_() {

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

    fprintf(stderr, "NODES = %zu, HITS = %zu\n", nodes_.size(), edges_.size());
}

Graph::~Graph() {
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

    uint32_t ttot = 0;
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
                    if (isSimilar(edge_xy->length + edge_yz->length, candidate_edge[z]->length, kTransitiveEdgeEps)) {
                        candidate_edge[z]->mark = true;
                        candidate_edge[z]->pair->mark = true;
                        ttot += 2;
                    }
                }
            }
        }

        for (const auto& edge: node_x->suffix_edges) {
            candidate_edge[edge->end_node->id] = nullptr;
        }
    }

    fprintf(stderr, "Removed %u transitive edges\n", ttot);
    remove_marked_edges();

    ttot = 0;
    for (const auto& it: edges_) {
        if (it == nullptr) continue;
        ++ttot;
    }
    fprintf(stderr, "%u remaining edges\n", ttot);
}

void Graph::remove_long_edges() {

    uint32_t total = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr) continue;

        for (const auto& edge1: node->suffix_edges) {
            for (const auto& edge2: node->suffix_edges) {
                if (edge1->id == edge2->id || edge1->mark == true || edge2->mark == true) continue;
                if (edge1->matching_bases() > kMinMatchingBasesRatio * edge2->matching_bases()) {
                    edge2->mark = true;
                    edge2->pair->mark = true;
                    ++total;
                }
            }
        }
    }

    fprintf(stderr, "Total long edges removed = %d\n", total);

    remove_marked_edges();
}

void Graph::remove_tips() {

    for (const auto& node: nodes_) {
        if (node == nullptr || !node->is_tip()) continue;

        uint32_t num_removed_edges = 0;

        for (const auto& edge: node->suffix_edges) {
            if (edge->end_node->in_degree() > 1) {
                // fprintf(stderr, "Removing %d\n", edge->begin_node->id);
                edge->mark = true;
                edge->pair->mark = true;
                ++num_removed_edges;
            }
        }

        if (num_removed_edges == node->suffix_edges.size()) {
            node->mark = true;
            node->pair->mark = true;
        }

        remove_marked_edges();
    }

    remove_isolated_nodes();
}

void Graph::remove_cycles() {

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

    do {
        cycles.clear();
        for (const auto& node: nodes_) {
            if (node == nullptr) continue;
            if (indexes[node->id] == -1) {
                strong_connect(node->id);
            }
        }

        fprintf(stderr, "Number of cycles %zu\n", cycles.size());

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
        }

        remove_marked_edges();

    } while (cycles.size() != 0);
}

void Graph::remove_bubbles() {

    std::vector<uint32_t> distance(nodes_.size(), 0);
    std::vector<uint32_t> visited(nodes_.size(), 0);
    uint32_t visited_length = 0;
    std::vector<int32_t> predecessor(nodes_.size(), -1);
    std::deque<uint32_t> queue;

    auto extract_path = [&](std::vector<uint32_t>& dst, int32_t source, int32_t sink) -> void {
        int32_t curr_id = sink;
        while (curr_id != source) {
            dst.push_back(curr_id);
            curr_id = predecessor[curr_id];
        }
        dst.push_back(source);
        std::reverse(dst.begin(), dst.end());
    };

    auto is_valid_bubble = [](const std::vector<uint32_t>& path, const std::vector<uint32_t>& other_path) -> bool {
        std::set<uint32_t> node_set;
        for (uint32_t i = 1; i < path.size() - 1; ++i) node_set.insert(path[i]);
        for (uint32_t i = 1; i < other_path.size() - 1; ++i) node_set.insert(other_path[i]);
        if (node_set.count(path[0]) != 0 || node_set.count(path[path.size()-1]) != 0) {
            return false;
        }
        if (path.size() + other_path.size() - 4 != node_set.size()) {
            return false;
        }
        return true;
    };

    auto path_params = [&](const std::vector<uint32_t>& path, uint32_t& num_reads, uint32_t& matching_bases, double& quality) -> void {
        num_reads = 0;
        for (const auto& it: path) num_reads += nodes_[it]->unitig_size;
        for (const auto& edge: nodes_[path[0]]->suffix_edges) {
            if (edge->end_node->id == path[1]) {
                quality = edge->quality;
                matching_bases = edge->matching_bases();
                break;
            }
        }
        return;
    };

    auto remove_path = [&](const std::vector<uint32_t>& path) -> void {
        for (uint32_t i = 0; i < path.size() - 1; ++i) {
            const auto& node = nodes_[path[i]];
            if (i != 0 && (node->in_degree() > 1 || node->out_degree() > 1)) {
                // fprintf(stderr, "Breaking at %d\n", path[i]);
                break;
            }
            for (const auto& edge: node->suffix_edges) {
                if (edge->end_node->id == path[i + 1]) {
                    edge->mark = true;
                    edge->pair->mark = true;
                    break;
                }
            }
        }
        return;
    };

    uint32_t bubbles_popped = 0;

    std::vector<uint32_t> bubble_candidates;
    locate_bubble_sources(bubble_candidates);
    for (const auto& id: bubble_candidates) {
        const auto& node = nodes_[id];

    // for (const auto& node: nodes_) {
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
            fprintf(stderr, "Source = %u, sink = %u, sink_predecesors = [%u, %u]\n", source, sink, predecessor[sink], sink_other_predecesor);

            std::vector<uint32_t> path, other_path;
            extract_path(path, source, sink);
            other_path.push_back(sink);
            extract_path(other_path, source, sink_other_predecesor);

            fprintf(stderr, "Path 1:");
            for (const auto& it: path) fprintf(stderr, " %d", it);
            fprintf(stderr, "\n");

            fprintf(stderr, "Path 2:");
            for (const auto& it: other_path) fprintf(stderr, " %d", it);
            fprintf(stderr, "\n");

            if (!is_valid_bubble(path, other_path)) {
                fprintf(stderr, "Not valid bubble!\n");
            } else {
                uint32_t path_reads = 0, path_matching_bases = 0;
                double path_quality = 0;
                path_params(path, path_reads, path_matching_bases, path_quality);

                uint32_t other_path_reads = 0, other_path_matching_bases = 0;
                double other_path_quality = 0;
                path_params(other_path, other_path_reads, other_path_matching_bases,
                    other_path_quality);

                fprintf(stderr, "Path 1 = (%d, %d, %g) | Path 2 = (%d, %d, %g) | Worse path is ",
                    path_reads, path_matching_bases, path_quality,
                    other_path_reads, other_path_matching_bases, other_path_quality);

                if (path_reads > other_path_reads || (path_reads == other_path_reads && path_matching_bases > other_path_matching_bases)) {
                    fprintf(stderr, "2\n");
                    remove_path(other_path);
                } else {
                    fprintf(stderr, "1\n");
                    remove_path(path);
                }
                remove_marked_edges();
                ++bubbles_popped;
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
    fprintf(stderr, "Popped %d bubbles\n", bubbles_popped);
}

void Graph::create_unitigs() {

    uint32_t node_id = nodes_.size();
    std::vector<bool> visited(nodes_.size(), false);
    std::vector<std::unique_ptr<Node>> new_nodes;

    for (const auto& node: nodes_) {
        if (node == nullptr || visited[node->id] || node->is_junction()) continue;

        auto bnode = node.get();
        while (!bnode->is_junction()) {
            visited[bnode->id] = true;
            visited[bnode->pair->id] = true;
            if (bnode->in_degree() == 0 || bnode->prefix_edges.front()->begin_node->is_junction()) {
                break;
            }
            bnode = bnode->prefix_edges.front()->begin_node;
        }

        auto enode = node.get();
        while (!enode->is_junction()) {
            visited[enode->id] = true;
            visited[enode->pair->id] = true;
            if (enode->out_degree() == 0 || enode->suffix_edges.front()->end_node->is_junction()) {
                break;
            }
            enode = enode->suffix_edges.front()->end_node;
        }

        if (bnode->id == enode->id) {
            continue;
        }

        Node* unitig = new Node(node_id++, bnode, enode); // normal
        Node* _unitig = new Node(node_id++, enode->pair, bnode->pair); // reverse complement

        unitig->pair = _unitig;
        _unitig->pair = unitig;

        new_nodes.push_back(std::unique_ptr<Node>(unitig));
        new_nodes.push_back(std::unique_ptr<Node>(_unitig));

        // fprintf(stderr, "Unitig: %d -> %d && %d -> %d\n", bnode->id, enode->id, enode->pair->id, bnode->pair->id);
    }

    for (auto& node: new_nodes) {
        nodes_.push_back(std::move(node));
    }

    remove_marked_edges();
    remove_isolated_nodes();
}

void Graph::print_contigs() const {

    uint32_t contig_id = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr ||node->id % 2 == 0 || node->unitig_size < kMinUnitigSize) continue;
        fprintf(stdout, ">%d\n%s\n", contig_id++, node->sequence.c_str());
    }
}

void Graph::locate_bubble_sources(std::vector<uint32_t>& dst) {

    // 0 - unmarked, 1 - in que, 2 - marked
    std::vector<uint8_t> marks(nodes_.size(), false);
    std::deque<uint32_t> node_que;

    for (const auto& node: nodes_) {
        if (node == nullptr || marks[node->id] == true) {
            continue;
        }

        node_que.push_back(node->id);
        while (!node_que.empty()) {
            const auto& curr_node = nodes_[node_que.front()];
            node_que.pop_front();
            //fprintf(stderr, "%d, %d\n", curr_node->id, marks[curr_node->id]);

            if (marks[curr_node->id] != 2) {
                if (curr_node->out_degree() > 1) {
                    dst.emplace_back(curr_node->id);
                }
                marks[curr_node->id] = 2;
                marks[curr_node->pair->id] = 2;

                for (const auto& edge: curr_node->prefix_edges) {
                    if (marks[edge->begin_node->id] == 0) {
                        marks[edge->begin_node->id] = 1;
                        node_que.push_back(edge->begin_node->id);
                    }
                }
                for (const auto& edge: curr_node->suffix_edges) {
                    if (marks[edge->end_node->id] == 0) {
                        marks[edge->end_node->id] = 1;
                        node_que.push_back(edge->end_node->id);
                    }
                }
            }
        }
    }

    /*for (const auto& it: dst) {
        fprintf(stderr, "%d ", it);
    }
    fprintf(stderr, "\n");*/
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

    for (const auto& node: nodes_) {
        if (node == nullptr) continue;
        delete_edges(node->prefix_edges);
        delete_edges(node->suffix_edges);
    }

    for (auto& edge: edges_) {
        if (edge == nullptr) continue;
        if (edge->mark == true) {
            edge.reset();
        }
    }
}

void Graph::print_dot() const {

    printf("digraph 1 {\n");
    printf("    overlap = scalexy\n");

    for (const auto& node: nodes_) {
        if (node == nullptr) continue;

        printf("    %d [label = \"%u [%u] {%d} U:%d\"", node->id, node->id, node->length(), node->read_id, node->unitig_size);
        if (node->id % 2 == 1) {
            printf(", style = filled, fillcolor = brown1]\n");
            printf("    %d -> %d [style = dotted, arrowhead = none]\n", node->id, node->id - 1);
        } else {
            printf("]\n");
        }
    }

    for (const auto& edge: edges_) {
        if (edge == nullptr) continue;
        printf("    %d -> %d [label = \"%d, %g\"]\n", edge->begin_node->id, edge->end_node->id, edge->length, edge->quality);
    }

    printf("}\n");
}

void Graph::print_csv() const {
    for (const auto& node: nodes_) {
        if (node == nullptr || node->id % 2 == 0) continue;
        printf("%u [%u] {%d} U:%d,%u [%u] {%d} U:%d,0,-\n",
            node->id, node->length(), node->read_id, node->unitig_size,
            node->pair->id, node->pair->length(), node->pair->read_id, node->pair->unitig_size);
    }
    for (const auto& edge: edges_) {
        if (edge == nullptr) continue;
        printf("%u [%u] {%d} U:%d,%u [%u] {%d} U:%d,1,%d %g\n",
            edge->begin_node->id, edge->begin_node->length(), edge->begin_node->read_id, edge->begin_node->unitig_size,
            edge->end_node->id, edge->end_node->length(), edge->end_node->read_id, edge->end_node->unitig_size,
            edge->length, edge->quality);
    }
}

}
