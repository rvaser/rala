/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <deque>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>

#include "sequence.hpp"
#include "overlap.hpp"
#include "pile.hpp"
#include "graph.hpp"

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "logger/logger.hpp"

namespace rala {

constexpr uint32_t kChunkSize = 1024 * 1024 * 1024; // ~1GB

bool comparable(double a, double b, double eps) {
    return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
        (b >= a * (1 - eps) && b <= a * (1 + eps));
}

template<typename T>
void shrinkToFit(std::vector<T>& src, uint64_t begin) {

    uint64_t i = begin;
    for (uint64_t j = begin; i < src.size(); ++i) {
        if (src[i] != nullptr) {
            continue;
        }

        j = std::max(j, i);
        while (j < src.size() && src[j] == nullptr) {
            ++j;
        }

        if (j >= src.size()) {
            break;
        } else if (i != j) {
            std::swap(src[i], src[j]);
        }
    }
    if (i < src.size()) {
        src.resize(i);
    }
}

class Graph::Node {
public:
    // Sequence encapsulation
    Node(uint64_t id, uint64_t sequence_id, const std::string& name,
        const std::string& data);
    // Unitig
    Node(uint64_t id, Node* begin_node, Node* end_node);
    Node(const Node&) = delete;
    const Node& operator=(const Node&) = delete;

    ~Node();

    bool is_rc() const {
        return id_ & 1;
    }

    uint32_t length() const {
        return data_.size();
    }

    uint32_t indegree() const {
        return prefix_edges_.size();
    }

    uint32_t outdegree() const {
        return suffix_edges_.size();
    }

    bool is_junction() const {
        return outdegree() > 1 || indegree() > 1;
    }

    bool is_tip() const {
        return outdegree() > 0 && indegree() == 0 && sequence_ids_.size() < 6;
    }

    uint64_t id_;
    std::string name_;
    std::string data_;
    std::vector<Edge*> prefix_edges_;
    std::vector<Edge*> suffix_edges_;
    std::vector<uint64_t> sequence_ids_;
    bool is_first_rc_;
    bool is_last_rc_;
    bool is_marked_;
    Node* pair_;
};

class Graph::Edge {
public:
    // Overlap encapsulatipn
    Edge(uint64_t id, Node* begin_node, Node* end_node, uint32_t length);
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    ~Edge();

    std::string label() const {
        return begin_node_->data_.substr(0, length_);
    }

    uint64_t id_;
    Node* begin_node_;
    Node* end_node_;
    uint32_t length_;
    double weight_;
    bool is_marked_;
    Edge* pair_;
};

Graph::Node::Node(uint64_t id, uint64_t sequence_id, const std::string& name,
    const std::string& data)
        : id_(id), name_(name), data_(data), prefix_edges_(), suffix_edges_(),
        sequence_ids_(1, sequence_id), is_first_rc_(id & 1), is_last_rc_(id & 1),
        is_marked_(false), pair_() {
}

Graph::Node::Node(uint64_t id, Node* begin_node, Node* end_node)
        : id_(id), name_(), data_(), prefix_edges_(), suffix_edges_(),
        sequence_ids_(), is_marked_(false), pair_() {

    if (begin_node == nullptr) {
        fprintf(stderr, "[rala::Graph::Node::Node] error: missing begin node!\n");
        exit(1);
    }
    if (end_node == nullptr) {
        fprintf(stderr, "[rala::Graph::Node::Node] error: missing end node!\n");
        exit(1);
    }

    is_first_rc_ = begin_node->is_first_rc_;

    auto node = begin_node;
    while (true) {
        auto edge = node->suffix_edges_[0];

        data_ += edge->label();
        sequence_ids_.insert(sequence_ids_.end(),
            node->sequence_ids_.begin(),
            node->sequence_ids_.end());
        is_last_rc_ = node->is_last_rc_;

        node = edge->end_node_;
        if (node == end_node) {
            break;
        }
    }

    if (begin_node != end_node) {
        data_ += end_node->data_;
        sequence_ids_.insert(sequence_ids_.end(),
            end_node->sequence_ids_.begin(),
            end_node->sequence_ids_.end());
        is_last_rc_ = end_node->is_last_rc_;
    }
}

Graph::Node::~Node() {
}

Graph::Edge::Edge(uint64_t id, Node* begin_node, Node* end_node, uint32_t length)
        : id_(id), begin_node_(begin_node), end_node_(end_node), length_(length),
        weight_(0), is_marked_(false), pair_() {
}

Graph::Edge::~Edge() {
}

std::unique_ptr<Graph> createGraph(const std::string& sequences_path,
    const std::string& overlaps_path, uint32_t num_threads) {

    std::unique_ptr<bioparser::Parser<Sequence>> sparser = nullptr;
    std::unique_ptr<bioparser::Parser<Overlap>> oparser = nullptr;

    auto is_suffix = [](const std::string& src, const std::string& suffix) -> bool {
        if (src.size() < suffix.size()) {
            return false;
        }
        return src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    if (is_suffix(sequences_path, ".fasta") || is_suffix(sequences_path, ".fa") ||
        is_suffix(sequences_path, ".fasta.gz") || is_suffix(sequences_path, ".fa.gz")) {
        sparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            sequences_path);
    } else if (is_suffix(sequences_path, ".fastq") || is_suffix(sequences_path, ".fq") ||
        is_suffix(sequences_path, ".fastq.gz") || is_suffix(sequences_path, ".fq.gz")) {
        sparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            sequences_path);
    } else {
        fprintf(stderr, "[rala::createGraph] error: "
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!\n",
            sequences_path.c_str());
        exit(1);
    }

    if (is_suffix(overlaps_path, ".mhap") || is_suffix(overlaps_path, ".mhap.gz")) {
        oparser = bioparser::createParser<bioparser::MhapParser, Overlap>(
            overlaps_path);
    } else if (is_suffix(overlaps_path, ".paf") || is_suffix(overlaps_path, ".paf.gz")) {
        oparser = bioparser::createParser<bioparser::PafParser, Overlap>(
            overlaps_path);
    } else {
        fprintf(stderr, "[rala::createGraph] error: "
            "file %s has unsupported format extension (valid extensions: "
            ".mhap, .mhap.gz, .paf, .paf.gz)!\n", overlaps_path.c_str());
        exit(1);
    }

    return std::unique_ptr<Graph>(new Graph(std::move(sparser), std::move(oparser),
        num_threads));
}

Graph::Graph(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
    std::unique_ptr<bioparser::Parser<Overlap>> oparser,
    uint32_t num_threads)
        : sparser_(std::move(sparser)), name_to_id_(), piles_(),
        coverage_median_(0), oparser_(std::move(oparser)), is_valid_overlap_(),
        thread_pool_(thread_pool::createThreadPool(num_threads)),
        nodes_(), edges_(), marked_edges_(), transitive_edges_(),
        logger_(logger::createLogger()) {
}

Graph::~Graph() {
    (*logger_).total("[rala::Graph::] total =");
}

void Graph::initialize() {

    (*logger_)();

    // create piles and sequence name hash
    uint64_t num_sequences = 0;
    sparser_->reset();
    while (true) {
        std::vector<std::unique_ptr<Sequence>> sequences;
        auto status = sparser_->parse_objects(sequences, kChunkSize);

        for (uint64_t i = 0; i < sequences.size(); ++i, ++num_sequences) {
            name_to_id_[sequences[i]->name()] = num_sequences;
            piles_.emplace_back(createPile(num_sequences,
                sequences[i]->data().size()));
        }

        if (!status) {
            break;
        }
    }

    (*logger_)("[rala::Graph::initialize] loaded sequences");
    (*logger_)();

    // update piles
    std::vector<std::unique_ptr<Overlap>> overlaps;
    uint64_t num_overlaps = 0;

    auto remove_duplicate_overlaps = [&](uint64_t begin, uint64_t end) -> void {
        for (uint64_t i = begin; i < end; ++i) {
            if (overlaps[i] == nullptr) {
                continue;
            }
            if (overlaps[i]->a_id() == overlaps[i]->b_id()) {
                // repeat
                // TODO: treat tandem and normal repeats differently?
                // if (overlaps[i]->a_begin() < overlaps[i]->b_end() &&
                //     overlaps[i]->b_begin() < overlaps[i]->a_end()) {
                // piles_[overlaps[i]->a_id()]->add_repetitive_region(
                //     std::min(overlaps[i]->a_begin(), overlaps[i]->b_begin()),
                //     std::max(overlaps[i]->a_end(), overlaps[i]->b_end()));
                is_valid_overlap_[num_overlaps + i] = false;
                //overlaps[i].reset();
                continue;
            }
            for (uint64_t j = i + 1; j < end; ++j) {
                if (overlaps[j] == nullptr) {
                    continue;
                }
                if (overlaps[i]->b_id() != overlaps[j]->b_id()) {
                    continue;
                }

                // TODO: process dual overlaps similar to self overlaps?
                if (overlaps[i]->length() > overlaps[j]->length()) {
                    is_valid_overlap_[num_overlaps + j] = false;
                } else {
                    is_valid_overlap_[num_overlaps + i] = false;
                    break;
                }
            }
        }
    };

    std::vector<std::vector<uint32_t>> overlap_bounds(piles_.size());

    auto store_overlap_bounds = [&](uint64_t begin, uint64_t end) -> void {
        for (uint64_t i = begin; i < end; ++i) {
            if (overlaps[i] == nullptr) {
                continue;
            }

            overlap_bounds[overlaps[i]->a_id()].emplace_back(
                (overlaps[i]->a_begin() + 15) << 1);
            overlap_bounds[overlaps[i]->a_id()].emplace_back(
                (overlaps[i]->a_end() - 15) << 1 | 1);
            overlap_bounds[overlaps[i]->b_id()].emplace_back(
                (overlaps[i]->b_begin() + 15) << 1);
            overlap_bounds[overlaps[i]->b_id()].emplace_back(
                (overlaps[i]->b_end() - 15) << 1 | 1);
        }
    };

    oparser_->reset();
    while (true) {
        uint64_t l = overlaps.size();
        auto status = oparser_->parse_objects(overlaps, kChunkSize);

        is_valid_overlap_.resize(is_valid_overlap_.size() + overlaps.size() - l, true);

        uint64_t c = 0;
        for (uint64_t i = l; i < overlaps.size(); ++i) {
            if (overlaps[i]->transmute(piles_, name_to_id_) == false) {
                overlaps[i].reset();
                is_valid_overlap_[num_overlaps + i] = false;
                continue;
            }

            while (overlaps[c] == nullptr) {
                ++c;
            }
            if (overlaps[c]->a_id() != overlaps[i]->a_id()) {
                remove_duplicate_overlaps(c, i);
                store_overlap_bounds(c, i);
                c = i;
            }
        }
        if (!status) {
            remove_duplicate_overlaps(c, overlaps.size());
            store_overlap_bounds(c, overlaps.size());
            c = overlaps.size();
        }
        num_overlaps += c;

        {
            std::vector<std::unique_ptr<Overlap>> tmp;
            for (uint64_t i = c; i < overlaps.size(); ++i) {
                tmp.emplace_back(std::move(overlaps[i]));
            }
            overlaps.swap(tmp);
        }

        std::vector<std::future<void>> thread_futures;
        for (const auto& it: piles_) {
            thread_futures.emplace_back(thread_pool_->submit_task(
                [&](uint64_t i) -> void {
                    piles_[i]->add_layers(overlap_bounds[i]);
                    std::vector<uint32_t>().swap(overlap_bounds[i]);
                }, it->id()));
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }

        if (!status) {
            break;
        }
    }

    (*logger_)("[rala::Graph::initialize] loaded overlaps");
    (*logger_)();

    std::vector<std::future<void>> thread_futures;
    for (const auto& it: piles_) {
        if (it == nullptr) {
            continue;
        }

        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint64_t i) -> void {
                if (piles_[i]->find_valid_region() == false) {
                    piles_[i].reset();
                } else {
                    piles_[i]->find_median();
                    piles_[i]->find_chimeric_hills();
                    piles_[i]->find_chimeric_pits();
                }
            }, it->id()));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    uint64_t num_prefiltered_sequences = 0;
    for (const auto& it: piles_) {
        if (it == nullptr) {
            ++num_prefiltered_sequences;
        }
    }

    (*logger_)("[rala::Graph::initialize] prefiltered sequences");

    if (num_prefiltered_sequences == num_sequences) {
        fprintf(stderr, "[rala::Graph::initialize] error: filtered all sequences!\n");
        exit(1);
    }

    fprintf(stderr, "[rala::Graph::initialize] number of prefiltered sequences = %lu\n",
        num_prefiltered_sequences);
}

void Graph::construct(const std::string& sensitive_overlaps_path) {

    if (!piles_.empty()) {
        fprintf(stderr, "[rala::Graph::construct] warning: "
            "object already constructed!\n");
        return;
    }

    initialize();

    (*logger_)();

    // store overlaps
    std::vector<std::unique_ptr<Overlap>> overlaps, internals;
    uint64_t num_overlaps = 0;

    oparser_->reset();
    while (true) {
        uint64_t l = overlaps.size();
        auto status = oparser_->parse_objects(overlaps, kChunkSize);

        for (uint64_t i = l; i < overlaps.size(); ++i) {
            auto& it = overlaps[i];
            if (is_valid_overlap_[num_overlaps + i - l] == false ||
                it->transmute(piles_, name_to_id_) == false ||
                it->trim(piles_) == false) {
                it.reset();
                continue;
            }

            if (piles_[it->a_id()]->has_chimeric_hill()) {
                piles_[it->a_id()]->check_chimeric_hills(it);
            }
            if (piles_[it->b_id()]->has_chimeric_hill()) {
                piles_[it->b_id()]->check_chimeric_hills(it);
            }

            switch (it->type(piles_)) {
                case OverlapType::kX:
                    internals.emplace_back(std::move(it));
                    //it.reset();
                    break;
                case OverlapType::kB:
                    if (piles_[it->b_id()]->has_chimeric_region() == false) {
                        piles_[it->a_id()].reset();
                        it.reset();
                    }
                    break;
                case OverlapType::kA:
                    if (piles_[it->a_id()]->has_chimeric_region() == false) {
                        piles_[it->b_id()].reset();
                        it.reset();
                    }
                    break;
                default:
                    break;
            }

            if (it == nullptr) {
                continue;
            }
        }
        num_overlaps += overlaps.size() - l;

        shrinkToFit(overlaps, l);

        if (!status) {
            // check if all non valid overlaps are deleted
            for (auto& it: overlaps) {
                if (it == nullptr) {
                    continue;
                }
                if (piles_[it->a_id()] == nullptr ||
                    piles_[it->b_id()] == nullptr) {

                    it.reset();
                }
            }
            shrinkToFit(overlaps, 0);

            for (auto& it: internals) {
                if (piles_[it->a_id()] == nullptr ||
                    piles_[it->b_id()] == nullptr) {

                    it.reset();
                }
            }
            shrinkToFit(internals, 0);

            break;
        }
    }

    (*logger_)("[rala::Graph::construct] loaded overlaps");

    preprocess(overlaps, internals);
    preprocess(overlaps, sensitive_overlaps_path);

    (*logger_)();

    // store reads
    std::vector<std::unique_ptr<Sequence>> sequences;

    sparser_->reset();
    while (true) {
        uint64_t l = sequences.size();
        auto status = sparser_->parse_objects(sequences, kChunkSize);

        for (uint64_t i = l; i < sequences.size(); ++i) {
            if (piles_[i] == nullptr) {
                sequences[i].reset();
                continue;
            }
            sequences[i]->trim(piles_[i]->begin(), piles_[i]->end());
            // piles_[i].reset();
        }

        if (!status) {
            break;
        }
    }

    (*logger_)("[rala::Graph::construct] loaded sequences");
    (*logger_)();

    // create assembly graph
    std::vector<int64_t> sequence_id_to_node_id(sequences.size(), -1);
    uint64_t node_id = 0;
    for (uint64_t i = 0; i < sequences.size(); ++i) {
        if (sequences[i] == nullptr) {
            continue;
        }
        const auto& it = sequences[i];

        sequence_id_to_node_id[i] = node_id;

        std::unique_ptr<Node> node(new Node(node_id++, i, it->name(), it->data()));
        std::unique_ptr<Node> node_complement(new Node(node_id++, i, it->name(),
            it->reverse_complement()));

        node->pair_ = node_complement.get();
        node_complement->pair_ = node.get();

        nodes_.emplace_back(std::move(node));
        nodes_.emplace_back(std::move(node_complement));

        sequences[i].reset();
    }

    uint64_t edge_id = 0;
    for (auto& it: overlaps) {
        Node* node_a = nodes_[sequence_id_to_node_id[it->a_id()]].get();
        Node* node_b = nodes_[sequence_id_to_node_id[it->b_id()] +
            it->orientation()].get();

        uint32_t a_length = piles_[it->a_id()]->end() - piles_[it->a_id()]->begin();
        uint32_t a_begin = it->a_begin() - piles_[it->a_id()]->begin();
        uint32_t a_end = it->a_end() - piles_[it->a_id()]->begin();

        uint32_t b_length = piles_[it->b_id()]->end() - piles_[it->b_id()]->begin();
        uint32_t b_begin = it->orientation() == 0 ?
            it->b_begin() - piles_[it->b_id()]->begin() :
            b_length - it->b_end() + piles_[it->b_id()]->begin();
        uint32_t b_end = it->orientation() == 0 ?
            it->b_end() - piles_[it->b_id()]->begin() :
            b_length - it->b_begin() + piles_[it->b_id()]->begin();

        if (it->type(piles_) == OverlapType::kAB) {
            std::unique_ptr<Edge> edge(new Edge(edge_id++, node_a, node_b,
                a_begin - b_begin));
            std::unique_ptr<Edge> edge_complement(new Edge(edge_id++,
                node_b->pair_, node_a->pair_, (b_length - b_end) -
                (a_length - a_end)));

            edge->pair_ = edge_complement.get();
            edge_complement->pair_ = edge.get();

            node_a->suffix_edges_.emplace_back(edge.get());
            node_a->pair_->prefix_edges_.emplace_back(edge_complement.get());
            node_b->prefix_edges_.emplace_back(edge.get());
            node_b->pair_->suffix_edges_.emplace_back(edge_complement.get());

            edges_.emplace_back(std::move(edge));
            edges_.emplace_back(std::move(edge_complement));

        } else if (it->type(piles_) == OverlapType::kBA) {
            std::unique_ptr<Edge> edge(new Edge(edge_id++, node_b, node_a,
                b_begin - a_begin));
            std::unique_ptr<Edge> edge_complement(new Edge(edge_id++,
                node_a->pair_, node_b->pair_, (a_length - a_end) -
                (b_length - b_end)));

            edge->pair_ = edge_complement.get();
            edge_complement->pair_ = edge.get();

            node_b->suffix_edges_.emplace_back(edge.get());
            node_b->pair_->prefix_edges_.emplace_back(edge_complement.get());
            node_a->prefix_edges_.emplace_back(edge.get());
            node_a->pair_->suffix_edges_.emplace_back(edge_complement.get());

            edges_.emplace_back(std::move(edge));
            edges_.emplace_back(std::move(edge_complement));
        }

        it.reset();
    }

    (*logger_)("[rala::Graph::construct] created assembly graph");

    fprintf(stderr, "[rala::Graph::construct] number of nodes = %zu\n",
        nodes_.size());
    fprintf(stderr, "[rala::Graph::construct] number of edges = %zu\n",
        edges_.size());
}

void Graph::simplify() {

    (*logger_)();

    uint32_t num_transitive_edges = remove_transitive_edges();

    uint32_t num_tips = 0;
    uint32_t num_bubbles = 0;

    while (true) {
        uint32_t num_changes = remove_tips();
        num_tips += num_changes;
        num_changes += num_changes;

        uint32_t num_changes_part = remove_bubbles();
        num_bubbles += num_changes_part;
        num_changes += num_changes_part;

        if (num_changes == 0) {
            break;
        }
    }

    shrink(42);
    uint32_t num_long_edges = 0;
    for (uint32_t i = 0; i < 5; ++i) {
        postprocess();
        num_long_edges += remove_long_edges();
        num_tips += remove_tips();
    }

    while (true) {
        uint32_t num_changes = remove_tips();
        num_tips += num_changes;
        num_changes += num_changes;

        uint32_t num_changes_part = remove_bubbles();
        num_bubbles += num_changes_part;
        num_changes += num_changes_part;

        if (num_changes == 0) {
            break;
        }
    }

    (*logger_)("[rala::Graph::simplify]");

    fprintf(stderr, "[rala::Graph::simplify] number of transitive edges = %u\n",
        num_transitive_edges);
    fprintf(stderr, "[rala::Graph::simplify] number of tips = %u\n",
        num_tips);
    fprintf(stderr, "[rala::Graph::simplify] number of bubbles = %u\n",
        num_bubbles);
    fprintf(stderr, "[rala::Graph::simplify] number of long edges = %u\n",
        num_long_edges);
}

void Graph::preprocess(std::vector<std::unique_ptr<Overlap>>& overlaps,
    std::vector<std::unique_ptr<Overlap>>& internals) {

    (*logger_)();

    std::vector<std::future<void>> thread_futures;
    for (const auto& it: piles_) {
        if (it == nullptr) {
            continue;
        }
        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint64_t i) -> void {
                if (piles_[i]->has_chimeric_hill() &&
                    piles_[i]->break_over_chimeric_hills() == false) {
                    piles_[i].reset();
                }
            }, it->id()));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    for (auto& it: overlaps) {
        if (it->trim(piles_) == false) {
            it.reset();
            continue;
        }
    }
    shrinkToFit(overlaps, 0);

    for (auto& it: internals) {
        if (it->trim(piles_) == false) {
            it.reset();
            continue;
        }
    }
    shrinkToFit(internals, 0);

    while (true) {

        std::vector<std::vector<uint64_t>> connections(piles_.size());
        for (const auto& it: overlaps) {
            connections[it->a_id()].emplace_back(it->b_id());
            connections[it->b_id()].emplace_back(it->a_id());
        }

        std::vector<std::vector<uint64_t>> components;
        std::vector<bool> is_visited(piles_.size(), false);
        for (uint64_t i = 0; i < connections.size(); ++i) {
            if (connections[i].empty() || is_visited[i]) {
                continue;
            }

            components.resize(components.size() + 1);

            std::deque<uint64_t> que = { i };
            while (!que.empty()) {
                uint64_t j = que.front();
                que.pop_front();

                if (is_visited[j]) {
                    continue;
                }
                is_visited[j] = true;
                components.back().emplace_back(j);

                for (const auto& it: connections[j]) {
                    que.emplace_back(it);
                }
                std::vector<uint64_t>().swap(connections[j]);
            }
        }
        std::vector<std::vector<uint64_t>>().swap(connections);
        std::vector<bool>().swap(is_visited);

        for (const auto& component: components) {

            std::vector<uint16_t> medians;
            for (const auto& it: component) {
                medians.emplace_back(piles_[it]->median());
            }
            std::nth_element(medians.begin(), medians.begin() + medians.size() / 2,
                medians.end());
            uint16_t component_median = medians[medians.size() / 2];

            for (const auto& it: component) {
                thread_futures.emplace_back(thread_pool_->submit_task(
                    [&](uint64_t i) -> void {
                        if (piles_[i]->break_over_chimeric_pits(component_median) == false) {
                            piles_[i].reset();
                        }
                    }, it));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }
            thread_futures.clear();
        }

        bool is_changed = false;

        for (auto& it: overlaps) {
            if (it->trim(piles_) == false) {
                it.reset();
                is_changed = true;
            }
        }
        shrinkToFit(overlaps, 0);

        for (auto& it: internals) {
            if (it->trim(piles_) == false) {
                it.reset();
                continue;
            }

            switch (it->type(piles_)) {
                case OverlapType::kAB:
                case OverlapType::kBA:
                    overlaps.emplace_back(std::move(it));
                    break;
                default:
                    break;
            }
        }
        shrinkToFit(internals, 0);

        if (is_changed == false) {
            break;
        }
    }

    for (auto& it: overlaps) {
        if (piles_[it->a_id()] == nullptr || piles_[it->b_id()] == nullptr) {
            it.reset();
            continue;
        }
        switch (it->type(piles_)) {
            case OverlapType::kA:
                piles_[it->b_id()].reset();
                it.reset();
                break;
            case OverlapType::kB:
                piles_[it->a_id()].reset();
                it.reset();
                break;
            default:
                break;
        }
    }
    for (auto& it: internals) {
        if (piles_[it->a_id()] == nullptr || piles_[it->b_id()] == nullptr) {
            it.reset();
            continue;
        }
        switch (it->type(piles_)) {
            case OverlapType::kA:
                piles_[it->b_id()].reset();
                it.reset();
                break;
            case OverlapType::kB:
                piles_[it->a_id()].reset();
                it.reset();
                break;
            default:
                break;
        }
    }
    shrinkToFit(internals, 0);

    for (auto& it: overlaps) {
        if (it == nullptr) {
            continue;
        }
        if (piles_[it->a_id()] == nullptr || piles_[it->b_id()] == nullptr) {
            it.reset();
        }
    }
    shrinkToFit(overlaps, 0);

    (*logger_)("[rala::Graph::preprocess]");
}

void Graph::preprocess(std::vector<std::unique_ptr<Overlap>>& overlaps,
    const std::string& path) {

    if (path.empty()) {
        return;
    }

    (*logger_)();

    std::unique_ptr<bioparser::Parser<Overlap>> oparser = nullptr;

    auto is_suffix = [](const std::string& src, const std::string& suffix) -> bool {
        if (src.size() < suffix.size()) {
            return false;
        }
        return src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    if (is_suffix(path, ".mhap") || is_suffix(path, ".mhap.gz")) {
        oparser = bioparser::createParser<bioparser::MhapParser, Overlap>(path);
    } else if (is_suffix(path, ".paf") || is_suffix(path, ".paf.gz")) {
        oparser = bioparser::createParser<bioparser::PafParser, Overlap>(path);
    } else {
        fprintf(stderr, "[rala::preprocess] error: "
            "file %s has unsupported format extension (valid extensions: "
            ".mhap, .mhap.gz, .paf, .paf.gz)!\n", path.c_str());
        exit(1);
    }

    uint64_t num_sequences = 0;
    std::unordered_map<uint64_t, uint64_t> sequence_id_to_id;
    std::vector<std::vector<uint32_t>> overlap_bounds;
    std::vector<std::unique_ptr<Overlap>> sensitive_overlaps;

    std::vector<std::future<void>> thread_futures;
    while (true) {
        uint64_t l = sensitive_overlaps.size();
        auto status = oparser->parse_objects(sensitive_overlaps, kChunkSize);

        for (uint64_t i = l; i < sensitive_overlaps.size(); ++i) {
            sensitive_overlaps[i]->transmute_(piles_, name_to_id_);
            if (sequence_id_to_id.find(sensitive_overlaps[i]->b_id()) == sequence_id_to_id.end()) {
                sequence_id_to_id[sensitive_overlaps[i]->b_id()] = num_sequences++;
            }
        }
        overlap_bounds.resize(num_sequences);

        for (uint64_t i = l; i < sensitive_overlaps.size(); ++i) {
            overlap_bounds[sequence_id_to_id[sensitive_overlaps[i]->b_id()]].emplace_back(
                sensitive_overlaps[i]->b_begin() << 1);
            overlap_bounds[sequence_id_to_id[sensitive_overlaps[i]->b_id()]].emplace_back(
                sensitive_overlaps[i]->b_end() << 1 | 1);

            if (sensitive_overlaps[i]->trim(piles_) == false) {
                sensitive_overlaps[i].reset();
            }
        }
        shrinkToFit(sensitive_overlaps, l);

        for (const auto& it: sequence_id_to_id) {
            if (!overlap_bounds[sequence_id_to_id[it.first]].empty()) {
                thread_futures.emplace_back(thread_pool_->submit_task(
                    [&](uint64_t i) -> void {
                        piles_[i]->add_layers(overlap_bounds[sequence_id_to_id[i]]);
                        std::vector<uint32_t>().swap(overlap_bounds[sequence_id_to_id[i]]);
                    }, it.first));
            }
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }
        thread_futures.clear();

        if (!status) {
            break;
        }
    }

    for (const auto& it: sequence_id_to_id) {
        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint64_t i) -> void {
                piles_[i]->find_median();
            }, it.first));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    std::vector<std::vector<uint64_t>> connections(piles_.size());
    for (const auto& it: overlaps) {
        connections[it->a_id()].emplace_back(it->b_id());
        connections[it->b_id()].emplace_back(it->a_id());
    }

    std::vector<std::vector<uint64_t>> components;
    std::vector<bool> is_visited(piles_.size(), false);
    for (uint64_t i = 0; i < connections.size(); ++i) {
        if (connections[i].empty() || is_visited[i]) {
            continue;
        }

        components.resize(components.size() + 1);

        std::deque<uint64_t> que = { i };
        while (!que.empty()) {
            uint64_t j = que.front();
            que.pop_front();

            if (is_visited[j]) {
                continue;
            }
            is_visited[j] = true;
            components.back().emplace_back(j);

            for (const auto& it: connections[j]) {
                que.emplace_back(it);
            }
            std::vector<uint64_t>().swap(connections[j]);
        }
    }
    std::vector<std::vector<uint64_t>>().swap(connections);
    std::vector<bool>().swap(is_visited);

    for (const auto& component: components) {

        std::vector<uint16_t> medians;
        for (const auto& it: component) {
            medians.emplace_back(piles_[it]->median());
        }
        std::nth_element(medians.begin(), medians.begin() + medians.size() / 2,
            medians.end());
        uint16_t component_median = medians[medians.size() / 2];

        for (const auto& it: component) {
            thread_futures.emplace_back(thread_pool_->submit_task(
                [&](uint64_t i) -> void {
                    piles_[i]->find_repetitive_hills(component_median);
                }, it));
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }
        thread_futures.clear();
    }

    for (auto& it: sensitive_overlaps) {
        if (it->trim(piles_) == false) {
            continue;
        }
        switch (it->type(piles_)) {
            case OverlapType::kAB:
            case OverlapType::kBA:
                if (piles_[it->b_id()]->has_repetitive_hills()) {
                    piles_[it->b_id()]->check_repetitive_hills(it);
                }
                break;
            default:
                it.reset();
                break;
        }
    }

    for (auto& it: overlaps) {
        if (!piles_[it->a_id()]->is_valid_overlap(it->a_begin(), it->a_end()) ||
            !piles_[it->b_id()]->is_valid_overlap(it->b_begin(), it->b_end())) {
            it.reset();
        }
    }
    shrinkToFit(overlaps, 0);

    (*logger_)("[rala::Graph::preprocess]");
}

void Graph::postprocess() {

    (*logger_)();

    if (transitive_edges_.empty() == false) {
        std::vector<std::pair<uint64_t, uint64_t>> tmp = { transitive_edges_[0] };
        for (uint64_t i = 1; i < transitive_edges_.size(); ++i) {
            if (nodes_[transitive_edges_[i].first] == nullptr ||
                nodes_[transitive_edges_[i].second] == nullptr) {
                continue;
            }
            if (transitive_edges_[i].first != transitive_edges_[i].second &&
                transitive_edges_[i] != transitive_edges_[i - 1]) {
                tmp.emplace_back(transitive_edges_[i]);
            }
        }
        tmp.swap(transitive_edges_);
    }

    std::vector<std::unordered_set<uint64_t>> components;
    std::vector<bool> is_visited(piles_.size(), false);
    for (uint64_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i] == nullptr || is_visited[i]) {
            continue;
        }

        components.resize(components.size() + 1);

        std::deque<uint64_t> que = { i };
        while (!que.empty()) {
            uint64_t j = que.front();
            que.pop_front();

            if (is_visited[j]) {
                continue;
            }
            const auto& node = nodes_[j];
            is_visited[node->id_] = true;
            is_visited[node->pair_->id_] = true;
            components.back().emplace((node->id_ >> 1) << 1);

            for (const auto& it: node->prefix_edges_) {
                que.emplace_back(it->begin_node_->id_);
            }
            for (const auto& it: node->suffix_edges_) {
                que.emplace_back(it->end_node_->id_);
            }
        }
    }
    std::vector<bool>().swap(is_visited);

    std::sort(components.begin(), components.end(),
        [](const std::unordered_set<uint64_t>& lhs, const std::unordered_set<uint64_t>& rhs) {
            return lhs.size() > rhs.size();
        }
    );

    std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<> distribution(0., 1.);

    using point = std::pair<double, double>;

    //uint32_t c = 0;
    for (const auto& component: components) {

        if (component.size() < 6) continue;

        bool has_junctions = false;
        for (const auto& it: component) {
            if (nodes_[it]->is_junction()) {
                has_junctions = true;
                break;
            }
        }
        if (has_junctions == false) continue;

        uint32_t num_iterations = 100;
        double k = sqrt(1. / static_cast<double>(component.size()));
        double t = 0.1;
        double dt = t / static_cast<double>(num_iterations + 1);

        auto add = [](const point& x, const point& y) {
            return std::make_pair(x.first + y.first, x.second + y.second);
        };
        auto substract = [](const point& x, const point& y) {
           return std::make_pair(x.first - y.first, x.second - y.second);
        };
        auto multiply = [](const point& x, double s) {
            return std::make_pair(x.first * s, x.second * s);
        };
        auto norm = [](const point& x) {
            return sqrt(x.first * x.first + x.second * x.second);
        };

        std::vector<point> points(nodes_.size());
        for (const auto& it: component) {
            points[it].first = distribution(generator);
            points[it].second = distribution(generator);
        }

        for (uint32_t i = 0; i < num_iterations; ++i) {
            std::vector<std::future<void>> thread_futures;
            std::vector<point> displacements(nodes_.size());

            auto thread_task = [&](uint64_t n) -> void {
                point displacement = {0., 0.};
                for (const auto& m: component) {
                    if (n == m) continue;
                    auto delta = substract(points[n], points[m]);
                    auto distance = norm(delta);
                    if (distance < 0.01) {
                        distance = 0.01;
                    }
                    displacement = add(displacement, multiply(delta, (k * k) / (distance * distance)));
                }
                for (const auto& e: nodes_[n]->prefix_edges_) {
                    auto m = (e->begin_node_->id_ >> 1) << 1;
                    auto delta = substract(points[n], points[m]);
                    auto distance = norm(delta);
                    if (distance < 0.01) {
                        distance = 0.01;
                    }
                    displacement = add(displacement, multiply(delta, -1. * distance / k));
                }
                for (const auto& e: nodes_[n]->suffix_edges_) {
                    auto m = (e->end_node_->id_ >> 1) << 1;
                    auto delta = substract(points[n], points[m]);
                    auto distance = norm(delta);
                    if (distance < 0.01) {
                        distance = 0.01;
                    }
                    displacement = add(displacement, multiply(delta, -1. * distance / k));
                }
                bool found = false;
                for (const auto& e: transitive_edges_) {
                    if (e.first != n) {
                        if (found) break;
                        continue;
                    }
                    found = true;
                    auto m = e.second;

                    auto delta = substract(points[n], points[m]);
                    auto distance = norm(delta);
                    if (distance < 0.01) {
                        distance = 0.01;
                    }
                    displacement = add(displacement, multiply(delta, -1. * distance / k));
                }

                auto length = norm(displacement);
                if (length < 0.01) {
                    length = 0.1;
                }
                displacements[n] = add(displacements[n], multiply(displacement, t / length));
                return;
            };

            for (const auto& n: component) {
                thread_futures.emplace_back(thread_pool_->submit_task(thread_task, n));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }
            for (const auto& n: component) {
                points[n] = add(points[n], displacements[n]);
            }

            t -= dt;
            ++i;
        }

        for (const auto& it: edges_) {
            if (it == nullptr || it->id_ & 1) {
                continue;
            }
            auto n = (it->begin_node_->id_ >> 1) << 1;
            auto m = (it->end_node_->id_ >> 1) << 1;

            if (component.find(n) != component.end() &&
                component.find(m) != component.end()) {
                it->weight_ = norm(substract(points[n], points[m]));
                it->pair_->weight_ = it->weight_;
            }
        }

        /*
        std::ofstream es("e" + std::to_string(c) + ".csv");
        std::ofstream ps("p" + std::to_string(c) + ".csv");
        std::ofstream ts("t" + std::to_string(c) + ".csv");
        ++c;

        for (const auto& it: component) {
            ps << it << "," << points[it].first << "," << points[it].second <<
                "," << (nodes_[it]->is_junction() ? 1 : 0) <<
                "," << nodes_[it]->sequence_ids_.size() << std::endl;
            for (const auto& e: nodes_[it]->prefix_edges_) {
                auto o = (e->begin_node_->id_ >> 1) << 1;
                es << it << "," << o << std::endl;
            }
            for (const auto& e: nodes_[it]->suffix_edges_) {
                auto o = (e->end_node_->id_ >> 1) << 1;
                es << it << "," << o << std::endl;
            }
        }

        for (const auto& e: transitive_edges_) {
            if (e.first < e.second &&
                component.find(e.first) != component.end() &&
                component.find(e.second) != component.end()) {
                ts << e.first << "," << e.second << std::endl;
            }
        }

        ts.close();
        ps.close();
        es.close();
        */
    }

    (*logger_)("[rala::Graph::postprocess]");

    return;
}

uint32_t Graph::remove_transitive_edges() {

    uint32_t num_transitive_edges = 0;
    std::vector<Edge*> candidate_edge(nodes_.size(), nullptr);

    for (const auto& node_a: nodes_) {
        if (node_a == nullptr) {
            continue;
        }

        for (const auto& edge_ab: node_a->suffix_edges_) {
            candidate_edge[edge_ab->end_node_->id_] = edge_ab;
        }

        for (const auto& edge_ab: node_a->suffix_edges_) {
            const auto& node_b = nodes_[edge_ab->end_node_->id_];

            for (const auto& edge_bc: node_b->suffix_edges_) {
                uint64_t c = edge_bc->end_node_->id_;

                if (candidate_edge[c] != nullptr && !candidate_edge[c]->is_marked_) {
                    if (comparable(edge_ab->length_ + edge_bc->length_,
                        candidate_edge[c]->length_, 0.12)) {

                        candidate_edge[c]->is_marked_ = true;
                        candidate_edge[c]->pair_->is_marked_ = true;
                        marked_edges_.emplace(candidate_edge[c]->id_);
                        marked_edges_.emplace(candidate_edge[c]->pair_->id_);
                        ++num_transitive_edges;
                    }
                }
            }
        }

        for (const auto& edge_ab: node_a->suffix_edges_) {
            candidate_edge[edge_ab->end_node_->id_] = nullptr;
        }
    }

    for (const auto& it: marked_edges_) {
        if (it & 1) {
            transitive_edges_.emplace_back(
                (edges_[it]->begin_node_->id_ >> 1) << 1,
                (edges_[it]->end_node_->id_ >> 1) << 1);
            transitive_edges_.emplace_back(
                transitive_edges_.back().second,
                transitive_edges_.back().first);
        }
    }
    std::sort(transitive_edges_.begin(), transitive_edges_.end());

    remove_marked_objects();

    return num_transitive_edges;
}

uint32_t Graph::remove_long_edges() {

    uint32_t num_long_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || node->suffix_edges_.size() < 2){
            continue;
        }

        for (const auto& edge: node->suffix_edges_) {
            for (const auto& other_edge: node->suffix_edges_) {
                if (edge->id_ == other_edge->id_ || edge->is_marked_ ||
                    other_edge->is_marked_) {
                    continue;
                }
                if (edge->weight_ * 2.0 < other_edge->weight_) {
                    other_edge->is_marked_ = true;
                    other_edge->pair_->is_marked_ = true;
                    marked_edges_.emplace(other_edge->id_);
                    marked_edges_.emplace(other_edge->pair_->id_);
                    ++num_long_edges;
                }
            }
        }
    }

    remove_marked_objects();

    return num_long_edges;
}

uint32_t Graph::remove_tips() {

    uint32_t num_tip_edges = 0;

    std::vector<bool> is_visited(nodes_.size(), false);

    for (const auto& it: nodes_) {
        if (it == nullptr || is_visited[it->id_] || !it->is_tip()) {
            continue;
        }

        bool is_circular = false;
        uint32_t num_reads = 0;

        auto end_node = it.get();
        while (!end_node->is_junction()) {
            num_reads += end_node->sequence_ids_.size();
            is_visited[end_node->id_] = true;
            is_visited[end_node->pair_->id_] = true;
            if (end_node->outdegree() == 0 ||
                end_node->suffix_edges_[0]->end_node_->is_junction()) {
                break;
            }
            end_node = end_node->suffix_edges_[0]->end_node_;
            if (end_node->id_ == it->id_) {
                is_circular = true;
                break;
            }
        }

        if (is_circular || end_node->outdegree() == 0 || num_reads > 5) {
            continue;
        }

        uint32_t num_removed_edges = 0;

        for (const auto& edge: end_node->suffix_edges_) {
            if (edge->end_node_->indegree() > 1) {
                edge->is_marked_ = true;
                edge->pair_->is_marked_ = true;
                marked_edges_.emplace(edge->id_);
                marked_edges_.emplace(edge->pair_->id_);
                ++num_removed_edges;
            }
        }

        if (num_removed_edges == end_node->suffix_edges_.size()) {
            auto curr_node = it.get();
            while (curr_node->id_ != end_node->id_) {
                curr_node->suffix_edges_[0]->is_marked_ = true;
                curr_node->suffix_edges_[0]->pair_->is_marked_ = true;
                marked_edges_.emplace(curr_node->suffix_edges_[0]->id_);
                marked_edges_.emplace(curr_node->suffix_edges_[0]->pair_->id_);
                curr_node = curr_node->suffix_edges_[0]->end_node_;
            }
        }

        num_tip_edges += num_removed_edges;

        remove_marked_objects(true);
    }

    return num_tip_edges;
}

uint32_t Graph::remove_bubbles() {

    std::vector<uint32_t> distance(nodes_.size(), 0);
    std::vector<uint64_t> visited(nodes_.size(), 0);
    uint64_t visited_length = 0;
    std::vector<int64_t> predecessor(nodes_.size(), -1);
    std::deque<uint64_t> node_queue;

    auto extract_path = [&](std::vector<uint64_t>& dst, uint64_t source,
        uint64_t sink) -> void {

        uint64_t curr_id = sink;
        while (curr_id != source) {
            dst.emplace_back(curr_id);
            curr_id = predecessor[curr_id];
        }
        dst.emplace_back(source);
        std::reverse(dst.begin(), dst.end());
    };

    auto calculate_path_length = [&](const std::vector<uint64_t>& path)
        -> uint32_t {

        if (path.empty()) {
            return 0;
        }

        uint32_t path_length = nodes_[path.back()]->length();
        for (uint64_t i = 0; i < path.size() - 1; ++i) {
            for (const auto& edge: nodes_[path[i]]->suffix_edges_) {
                if (edge->end_node_->id_ == (uint64_t) path[i + 1]) {
                    path_length += edge->length_;
                    break;
                }
            }
        }
        return path_length;
    };

    auto is_valid_bubble = [&](const std::vector<uint64_t>& path,
        const std::vector<uint64_t>& other_path) -> bool {

        if (path.empty() || other_path.empty()) {
            return false;
        }

        std::unordered_set<uint64_t> node_set;
        for (const auto& it: path) {
            node_set.emplace(it);
        }
        for (const auto& it: other_path) {
            node_set.emplace(it);
        }
        if (path.size() + other_path.size() - 2 != node_set.size()) {
            return false;
        }
        for (const auto& it: path) {
            uint64_t pair_id = nodes_[it]->pair_->id_;
            if (node_set.count(pair_id) != 0) {
                return false;
            }
        }
        uint32_t path_length = calculate_path_length(path);
        uint32_t other_path_length = calculate_path_length(other_path);
        if (std::min(path_length, other_path_length) <
            std::max(path_length, other_path_length) * 0.8) {

            for (uint64_t i = 1; i < other_path.size() - 1; ++i) {
                if (nodes_[other_path[i]]->is_junction()) {
                    return false;
                }
            }
            for (uint64_t i = 1; i < path.size() - 1; ++i) {
                if (nodes_[path[i]]->is_junction()) {
                    return false;
                }
            }
        }
        return true;
    };

    uint32_t num_bubbles_popped = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr || node->outdegree() < 2) {
            continue;
        }

        bool found_sink = false;
        uint64_t sink = 0, sink_other_predecesor = 0;
        uint64_t source = node->id_;

        // BFS
        node_queue.emplace_back(source);
        visited[visited_length++] = source;
        while (!node_queue.empty() && !found_sink) {
            uint64_t v = node_queue.front();
            const auto& curr_node = nodes_[v];

            node_queue.pop_front();

            for (const auto& edge: curr_node->suffix_edges_) {
                uint64_t w = edge->end_node_->id_;

                if (w == source) {
                    // Cycle
                    continue;
                }

                if (distance[v] + edge->length_ > 5000000) {
                    // Out of reach
                    continue;
                }

                distance[w] = distance[v] + edge->length_;
                visited[visited_length++] = w;
                node_queue.emplace_back(w);

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
            std::vector<uint64_t> path;
            extract_path(path, source, sink);

            std::vector<uint64_t> other_path(1, sink);
            extract_path(other_path, source, sink_other_predecesor);

            if (is_valid_bubble(path, other_path)) {
                uint64_t path_num_reads = 0;
                for (const auto& it: path) {
                    path_num_reads += nodes_[it]->sequence_ids_.size();
                }

                uint64_t other_path_num_reads = 0;
                for (const auto& it: other_path) {
                    other_path_num_reads += nodes_[it]->sequence_ids_.size();
                }

                std::vector<uint64_t> edges_for_removal;
                if (path_num_reads > other_path_num_reads) {
                    find_removable_edges(edges_for_removal, other_path);
                } else {
                    find_removable_edges(edges_for_removal, path);
                }

                if (edges_for_removal.empty()) {
                    uint32_t path_length = calculate_path_length(path);
                    uint32_t other_path_length = calculate_path_length(other_path);
                    if (std::min(path_length, other_path_length) >=
                        std::max(path_length, other_path_length) * 0.8) {

                        if (path_num_reads > other_path_num_reads) {
                            find_removable_edges(edges_for_removal, path);
                        } else {
                            find_removable_edges(edges_for_removal, other_path);
                        }
                    }
                }

                for (const auto& edge_id: edges_for_removal) {
                    edges_[edge_id]->is_marked_ = true;
                    edges_[edge_id]->pair_->is_marked_ = true;
                    marked_edges_.emplace(edge_id);
                    marked_edges_.emplace(edges_[edge_id]->pair_->id_);
                }
                if (!edges_for_removal.empty()) {
                    remove_marked_objects(true);
                    ++num_bubbles_popped;
                }
            }
        }

        node_queue.clear();
        for (uint64_t i = 0; i < visited_length; ++i) {
            distance[visited[i]] = 0;
            predecessor[visited[i]] = -1;
        }
        visited_length = 0;
    }

    return num_bubbles_popped;
}

uint64_t Graph::find_edge(uint64_t src, uint64_t dst) {

    uint64_t edge_id = 0;
    bool found_edge = false;
    for (const auto& edge: nodes_[src]->suffix_edges_) {
        if (edge->end_node_->id_ == dst) {
            edge_id = edge->id_;
            found_edge = true;
            break;
        }
    }

    if (!found_edge) {
        fprintf(stderr, "[rala::Graph::find_edge] error: "
            "missing edge between nodes %lu and %lu\n", src, dst);
        exit(1);
    }

    return edge_id;
}

void Graph::find_removable_edges(std::vector<uint64_t>& dst,
    const std::vector<uint64_t>& path) {

    if (path.empty()) {
        return;
    }

    // find first node with multiple in edges
    int64_t pref = -1;
    for (uint64_t i = 1; i < path.size() - 1; ++i) {
        if (nodes_[path[i]]->indegree() > 1) {
            pref = i;
            break;
        }
    }
    // find last node with multiple out edges
    int64_t suff = -1;
    for (uint64_t i = 1; i < path.size() - 1; ++i) {
        if (nodes_[path[i]]->outdegree() > 1) {
            suff = i;
        }
    }

    if (pref == -1 && suff == -1) {
        // remove whole path
        for (uint64_t i = 0; i < path.size() - 1; ++i) {
            dst.emplace_back(find_edge(path[i], path[i + 1]));
        }
        return;
    }

    if (pref != -1 && nodes_[path[pref]]->outdegree() > 1) {
        return;
    }
    if (suff != -1 && nodes_[path[suff]]->indegree() > 1) {
        return;
    }

    if (pref == -1) {
        // remove everything after last suff node
        for (uint64_t i = suff; i < path.size() - 1; ++i) {
            dst.emplace_back(find_edge(path[i], path[i + 1]));
        }
    } else if (suff == -1) {
        // remove everything before first pref node
        for (int64_t i = 0; i < pref; ++i) {
            dst.emplace_back(find_edge(path[i], path[i + 1]));
        }
    } else if (suff < pref) {
        // remove everything between last suff and first pref node
        for (int64_t i = suff; i < pref; ++i) {
            dst.emplace_back(find_edge(path[i], path[i + 1]));
        }
    }
}

uint32_t Graph::create_unitigs() {

    std::vector<bool> is_visited(nodes_.size(), false);

    uint64_t node_id = nodes_.size();
    std::vector<std::unique_ptr<Node>> unitigs;

    uint64_t edge_id = edges_.size();
    std::vector<std::unique_ptr<Edge>> unitig_edges;

    uint32_t num_unitigs_created = 0;

    for (const auto& it: nodes_) {
        if (it == nullptr || is_visited[it->id_] || it->is_junction()) {
            continue;
        }

        bool is_circular = false;
        auto begin_node = it.get();
        while (!begin_node->is_junction()) {
            is_visited[begin_node->id_] = true;
            is_visited[begin_node->pair_->id_] = true;
            if (begin_node->indegree() == 0 ||
                begin_node->prefix_edges_[0]->begin_node_->is_junction()) {
                break;
            }
            begin_node = begin_node->prefix_edges_[0]->begin_node_;
            if (begin_node->id_ == it->id_) {
                is_circular = true;
                break;
            }
        }

        auto end_node = it.get();
        while (!end_node->is_junction()) {
            is_visited[end_node->id_] = true;
            is_visited[end_node->pair_->id_] = true;
            if (end_node->outdegree() == 0 ||
                end_node->suffix_edges_[0]->end_node_->is_junction()) {
                break;
            }
            end_node = end_node->suffix_edges_[0]->end_node_;
            if (end_node->id_ == it->id_) {
                is_circular = true;
                break;
            }
        }

        if (!is_circular && begin_node == end_node) {
            continue;
        }

        std::unique_ptr<Node> unitig(new Node(node_id++, begin_node, end_node));
        std::unique_ptr<Node> unitig_complement(new Node(node_id++,
            end_node->pair_, begin_node->pair_));

        unitig->pair_ = unitig_complement.get();
        unitig_complement->pair_ = unitig.get();

        if (begin_node != end_node) {
            if (begin_node->indegree() != 0) {
                const auto& edge = begin_node->prefix_edges_[0];

                edge->is_marked_ = true;
                edge->pair_->is_marked_ = true;
                marked_edges_.emplace(edge->id_);
                marked_edges_.emplace(edge->pair_->id_);

                std::unique_ptr<Edge> unitig_edge(new Edge(edge_id++,
                    edge->begin_node_, unitig.get(), edge->length_));
                std::unique_ptr<Edge> unitig_edge_complement(new Edge(edge_id++,
                    unitig_complement.get(), edge->pair_->end_node_,
                    edge->pair_->length_ + unitig_complement->length() -
                    begin_node->pair_->length()));

                unitig_edge->pair_ = unitig_edge_complement.get();
                unitig_edge_complement->pair_ = unitig_edge.get();

                edge->begin_node_->suffix_edges_.emplace_back(unitig_edge.get());
                edge->pair_->end_node_->prefix_edges_.emplace_back(
                    unitig_edge_complement.get());
                unitig->prefix_edges_.emplace_back(unitig_edge.get());
                unitig_complement->suffix_edges_.emplace_back(
                    unitig_edge_complement.get());

                unitig_edges.emplace_back(std::move(unitig_edge));
                unitig_edges.emplace_back(std::move(unitig_edge_complement));
            }

            if (end_node->outdegree() != 0) {
                const auto& edge = end_node->suffix_edges_[0];

                edge->is_marked_ = true;
                edge->pair_->is_marked_ = true;
                marked_edges_.emplace(edge->id_);
                marked_edges_.emplace(edge->pair_->id_);

                std::unique_ptr<Edge> unitig_edge(new Edge(edge_id++,
                    unitig.get(), edge->end_node_, edge->length_ +
                    unitig->length() - end_node->length()));
                std::unique_ptr<Edge> unitig_edge_complement(new Edge(edge_id++,
                    edge->pair_->begin_node_, unitig_complement.get(),
                    edge->pair_->length_));

                unitig_edge->pair_ = unitig_edge_complement.get();
                unitig_edge_complement->pair_ = unitig_edge.get();

                unitig->suffix_edges_.emplace_back(unitig_edge.get());
                unitig_complement->prefix_edges_.emplace_back(
                    unitig_edge_complement.get());
                edge->end_node_->prefix_edges_.emplace_back(unitig_edge.get());
                edge->pair_->begin_node_->suffix_edges_.emplace_back(
                    unitig_edge_complement.get());

                unitig_edges.emplace_back(std::move(unitig_edge));
                unitig_edges.emplace_back(std::move(unitig_edge_complement));
            }
        }

        unitigs.emplace_back(std::move(unitig));
        unitigs.emplace_back(std::move(unitig_complement));

        ++num_unitigs_created;

        // mark edges for deletion
        auto node = begin_node;
        while (true) {
            const auto& edge = node->suffix_edges_[0];

            edge->is_marked_ = true;
            edge->pair_->is_marked_ = true;
            marked_edges_.emplace(edge->id_);
            marked_edges_.emplace(edge->pair_->id_);

            node = edge->end_node_;
            if (node == end_node) {
                break;
            }
        }
    }

    for (uint64_t i = 0; i < unitigs.size(); ++i) {
        nodes_.emplace_back(std::move(unitigs[i]));
    }
    for (uint64_t i = 0; i < unitig_edges.size(); ++i) {
        edges_.emplace_back(std::move(unitig_edges[i]));
    }

    remove_marked_objects(true);

    return num_unitigs_created;
}

uint32_t Graph::shrink(uint32_t epsilon) {

    std::vector<bool> is_visited(nodes_.size(), false);
    std::vector<uint64_t> node_updates(nodes_.size(), 0);

    uint64_t node_id = nodes_.size();
    std::vector<std::unique_ptr<Node>> unitigs;

    uint64_t edge_id = edges_.size();
    std::vector<std::unique_ptr<Edge>> unitig_edges;

    uint32_t num_unitigs_created = 0;

    for (const auto& it: nodes_) {
        if (it == nullptr || is_visited[it->id_] || it->is_junction()) {
            continue;
        }

        uint32_t extension = 1;

        bool is_circular = false;
        auto begin_node = it.get();
        while (!begin_node->is_junction()) {
            is_visited[begin_node->id_] = true;
            is_visited[begin_node->pair_->id_] = true;
            if (begin_node->indegree() == 0 ||
                begin_node->prefix_edges_[0]->begin_node_->is_junction()) {
                break;
            }
            begin_node = begin_node->prefix_edges_[0]->begin_node_;
            ++extension;
            if (begin_node->id_ == it->id_) {
                is_circular = true;
                break;
            }
        }

        if (is_circular) {
            continue;
        }

        auto end_node = it.get();
        while (!end_node->is_junction()) {
            is_visited[end_node->id_] = true;
            is_visited[end_node->pair_->id_] = true;
            if (end_node->outdegree() == 0 ||
                end_node->suffix_edges_[0]->end_node_->is_junction()) {
                break;
            }
            end_node = end_node->suffix_edges_[0]->end_node_;
            ++extension;
            if (end_node->id_ == it->id_) {
                is_circular = true;
                break;
            }
        }

        if (is_circular || begin_node == end_node || extension < 2 * epsilon + 2) {
            continue;
        }

        // update begin_node
        for (uint32_t i = 0; i < epsilon; ++i) {
            begin_node = begin_node->suffix_edges_[0]->end_node_;
        }

        // update end_node
        for (uint32_t i = 0; i < epsilon; ++i) {
            end_node = end_node->prefix_edges_[0]->begin_node_;
        }

        // update node ids for transitive edges
        auto node = begin_node;
        while (node != end_node) {
            node_updates[(node->id_ >> 1) << 1] = node_id;
            node = node->suffix_edges_[0]->end_node_;
        }

        std::unique_ptr<Node> unitig(new Node(node_id++, begin_node, end_node));
        std::unique_ptr<Node> unitig_complement(new Node(node_id++,
            end_node->pair_, begin_node->pair_));

        unitig->pair_ = unitig_complement.get();
        unitig_complement->pair_ = unitig.get();

        if (begin_node->indegree() != 0) {
            const auto& edge = begin_node->prefix_edges_[0];

            edge->is_marked_ = true;
            edge->pair_->is_marked_ = true;
            marked_edges_.emplace(edge->id_);
            marked_edges_.emplace(edge->pair_->id_);

            std::unique_ptr<Edge> unitig_edge(new Edge(edge_id++,
                edge->begin_node_, unitig.get(), edge->length_));
            std::unique_ptr<Edge> unitig_edge_complement(new Edge(edge_id++,
                unitig_complement.get(), edge->pair_->end_node_,
                edge->pair_->length_ + unitig_complement->length() -
                begin_node->pair_->length()));

            unitig_edge->pair_ = unitig_edge_complement.get();
            unitig_edge_complement->pair_ = unitig_edge.get();

            edge->begin_node_->suffix_edges_.emplace_back(unitig_edge.get());
            edge->pair_->end_node_->prefix_edges_.emplace_back(
                unitig_edge_complement.get());
            unitig->prefix_edges_.emplace_back(unitig_edge.get());
            unitig_complement->suffix_edges_.emplace_back(
                unitig_edge_complement.get());

            unitig_edges.emplace_back(std::move(unitig_edge));
            unitig_edges.emplace_back(std::move(unitig_edge_complement));
        }

        if (end_node->outdegree() != 0) {
            const auto& edge = end_node->suffix_edges_[0];

            edge->is_marked_ = true;
            edge->pair_->is_marked_ = true;
            marked_edges_.emplace(edge->id_);
            marked_edges_.emplace(edge->pair_->id_);

            std::unique_ptr<Edge> unitig_edge(new Edge(edge_id++,
                unitig.get(), edge->end_node_, edge->length_ +
                unitig->length() - end_node->length()));
            std::unique_ptr<Edge> unitig_edge_complement(new Edge(edge_id++,
                edge->pair_->begin_node_, unitig_complement.get(),
                edge->pair_->length_));

            unitig_edge->pair_ = unitig_edge_complement.get();
            unitig_edge_complement->pair_ = unitig_edge.get();

            unitig->suffix_edges_.emplace_back(unitig_edge.get());
            unitig_complement->prefix_edges_.emplace_back(
                unitig_edge_complement.get());
            edge->end_node_->prefix_edges_.emplace_back(unitig_edge.get());
            edge->pair_->begin_node_->suffix_edges_.emplace_back(
                unitig_edge_complement.get());

            unitig_edges.emplace_back(std::move(unitig_edge));
            unitig_edges.emplace_back(std::move(unitig_edge_complement));
        }

        unitigs.emplace_back(std::move(unitig));
        unitigs.emplace_back(std::move(unitig_complement));

        ++num_unitigs_created;

        // mark edges for deletion
        node = begin_node;
        while (true) {
            const auto& edge = node->suffix_edges_[0];

            edge->is_marked_ = true;
            edge->pair_->is_marked_ = true;
            marked_edges_.emplace(edge->id_);
            marked_edges_.emplace(edge->pair_->id_);

            node = edge->end_node_;
            if (node == end_node) {
                break;
            }
        }
    }

    for (uint64_t i = 0; i < unitigs.size(); ++i) {
        nodes_.emplace_back(std::move(unitigs[i]));
    }
    for (uint64_t i = 0; i < unitig_edges.size(); ++i) {
        edges_.emplace_back(std::move(unitig_edges[i]));
    }

    remove_marked_objects(true);

    // update transitive edges
    for (auto& it: transitive_edges_) {
        if (node_updates[it.first] != 0) {
            it.first = node_updates[it.first];
        }
        if (node_updates[it.second] != 0) {
            it.second = node_updates[it.second];
        }
    }
    std::sort(transitive_edges_.begin(), transitive_edges_.end());

    return num_unitigs_created;
}

void Graph::extract_contigs(std::vector<std::unique_ptr<Sequence>>& dst,
    bool drop_unassembled_sequences) {

    create_unitigs();

    uint32_t contig_id = 0;
    std::vector<uint32_t> contig_length;
    for (const auto& node: nodes_) {
        if (node == nullptr || node->is_rc()) {
            continue;
        }
        if (drop_unassembled_sequences && (node->sequence_ids_.size() < 6 ||
            node->length() < 10000)) {
            continue;
        }
        contig_length.emplace_back(node->data_.size());

        std::string name = "Ctg" + std::to_string(contig_id);
        name += " RC:i:" + std::to_string(node->sequence_ids_.size());
        name += " LN:i:" + std::to_string(node->data_.size());

        dst.emplace_back(createSequence(name, node->data_));
        ++contig_id;
    }

    fprintf(stderr, "[rala::Graph::extract_contigs] number of contigs = %zu\n",
        contig_length.size());

    if (contig_length.empty()) {
        return;
    }

    std::sort(contig_length.begin(), contig_length.end());

    fprintf(stderr, "[rala::Graph::extract_contigs] shortest contig length = %u\n",
        contig_length.front());
    fprintf(stderr, "[rala::Graph::extract_contigs] median contig length = %u\n",
        contig_length[contig_length.size() / 2]);
    fprintf(stderr, "[rala::Graph::extract_contigs] longest contig length = %u\n",
        contig_length.back());
}

void Graph::extract_nodes(std::vector<std::unique_ptr<Sequence>>& dst) {

    std::unordered_set<uint64_t> node_ids;
    for (const auto& it: nodes_) {
        if (it == nullptr || it->is_rc() || (it->outdegree() == 0 && it->indegree() == 0)) {
            continue;
        }

        node_ids.emplace(it->id_);

        for (const auto& edge: it->prefix_edges_) {
            if (edge->begin_node_->is_rc()) {
                node_ids.emplace(edge->begin_node_->pair_->id_);
            } else {
                node_ids.emplace(edge->begin_node_->id_);
            }
        }
        for (const auto& edge: it->suffix_edges_) {
            if (edge->end_node_->is_rc()) {
                node_ids.emplace(edge->end_node_->pair_->id_);
            } else {
                node_ids.emplace(edge->end_node_->id_);
            }
        }
    }

    for (const auto& it: node_ids) {
        dst.emplace_back(createSequence(nodes_[it]->name_, nodes_[it]->data_));
    }

    fprintf(stderr, "[rala::Graph::extract_nodes] number of nodes = %zu\n",
        dst.size());
}

void Graph::remove_marked_objects(bool remove_nodes) {

    auto delete_edges = [&](std::vector<Edge*>& edges) -> void {
        for (uint32_t i = 0; i < edges.size(); ++i) {
            if (edges[i]->is_marked_) {
                edges[i] = nullptr;
            }
        }
        shrinkToFit(edges, 0);
    };

    std::unordered_set<uint32_t> marked_nodes;
    for (const auto& it: marked_edges_) {
        if (remove_nodes) {
            marked_nodes.emplace(edges_[it]->begin_node_->id_);
            marked_nodes.emplace(edges_[it]->end_node_->id_);
        }
        delete_edges(edges_[it]->begin_node_->suffix_edges_);
        delete_edges(edges_[it]->end_node_->prefix_edges_);
    }

    if (remove_nodes) {
        for (const auto& it: marked_nodes) {
            if (nodes_[it]->outdegree() == 0 && nodes_[it]->indegree() == 0) {
                nodes_[it].reset();
            }
        }
    }

    for (const auto& it: marked_edges_) {
        edges_[it].reset();
    }
    marked_edges_.clear();
}

void Graph::print_csv(const std::string& path) const {

    auto graph_file = fopen(path.c_str(), "w");

    for (const auto& it: nodes_) {
        if (it == nullptr || !it->is_rc() || (it->outdegree() == 0 && it->indegree() == 0)) {
            continue;
        }
        fprintf(graph_file, "%lu LN:i:%u RC:i:%lu,%lu LN:i:%u RC:i:%lu,0,-\n",
            it->id_, it->length(), it->sequence_ids_.size(),
            it->pair_->id_, it->pair_->length(), it->pair_->sequence_ids_.size());
    }

    for (const auto& it: edges_) {
        if (it == nullptr) {
            continue;
        }
        fprintf(graph_file, "%lu LN:i:%u RC:i:%lu,%lu LN:i:%u RC:i:%lu,1,%lu %u %lf\n",
            it->begin_node_->id_, it->begin_node_->length(), it->begin_node_->sequence_ids_.size(),
            it->end_node_->id_, it->end_node_->length(), it->end_node_->sequence_ids_.size(),
            it->id_, it->length_, it->weight_);
    }

    fclose(graph_file);
}

void Graph::print_gfa(const std::string& path) const {

    auto graph_file = fopen(path.c_str(), "w");

    std::unordered_map<uint64_t, std::string> node_id_to_unitig_name;
    uint32_t unitig_id = 0;

    for (const auto& it: nodes_) {
        if (it == nullptr || it->is_rc() || (it->outdegree() == 0 && it->indegree() == 0)) {
            continue;
        }
        if (it->name_.empty()) {
            std::string unitig_name = "Utg" + std::to_string(unitig_id++);
            node_id_to_unitig_name[it->id_] = unitig_name;
            node_id_to_unitig_name[it->pair_->id_] = unitig_name;
        }

        const auto& node_name = !it->name_.empty() ? it->name_ :
            node_id_to_unitig_name[it->id_];

        fprintf(graph_file, "S\t%s\t%s\tLN:i:%zu\tRC:i:%lu\n",
            node_name.c_str(), it->data_.c_str(), it->data_.size(),
            it->sequence_ids_.size());
    }

    for (const auto& it: edges_) {
        if (it == nullptr) {
            continue;
        }

        const auto& begin_node_name = !nodes_[it->begin_node_->id_]->name_.empty() ?
            nodes_[it->begin_node_->id_]->name_ :
            node_id_to_unitig_name[it->begin_node_->id_];

        const auto& end_node_name = !nodes_[it->end_node_->id_]->name_.empty() ?
            nodes_[it->end_node_->id_]->name_ :
            node_id_to_unitig_name[it->end_node_->id_];

        fprintf(graph_file, "L\t%s\t%c\t%s\t%c\t%zuM\n",
            begin_node_name.c_str(), it->begin_node_->is_rc() ? '-' : '+',
            end_node_name.c_str(), it->end_node_->is_rc() ? '-' : '+',
            it->begin_node_->data_.size() - it->length_);
    }

    fclose(graph_file);
}

void Graph::print_json(const std::string& path) const {

    std::ofstream os(path);
    os << "{\"nodes\":{";
    bool is_first = true;

    std::unordered_set<uint64_t> sequence_ids;
    for (const auto& it: nodes_) {
        if (it == nullptr || it->is_rc() || !it->is_junction()) {
            continue;
        }

        if (!is_first) {
            os << ",";
        }
        is_first = false;
        os << "\"" << it->sequence_ids_.front() << "\":{\"n\":" << it->id_ << ",";

        os << "\"p\":[";

        sequence_ids.emplace(it->sequence_ids_.front());
        for (uint32_t i = 0; i < it->prefix_edges_.size(); ++i) {
            auto other = it->prefix_edges_[i]->begin_node_;
            sequence_ids.emplace(other->sequence_ids_.back());

            os << "[\"" << other->sequence_ids_.back() << "\",\"" <<
                other->id_ << "\"," << other->is_last_rc_ << "," <<
                other->length() - it->prefix_edges_[i]->length_ << "]";
            if (i < it->prefix_edges_.size() - 1) {
                os << ",";
            }
        }

        os << "],\"s\":[";

        for (uint32_t i = 0; i < it->suffix_edges_.size(); ++i) {
            auto other = it->suffix_edges_[i]->end_node_;
            sequence_ids.emplace(other->sequence_ids_.front());

            os << "[\"" << other->sequence_ids_.front() << "\",\"" <<
                other->id_ << "\"," << other->is_first_rc_ << "," <<
                it->length() - it->suffix_edges_[i]->length_ << "]";
            if (i < it->suffix_edges_.size() - 1) {
                os << ",";
            }
        }

        os << "]}";
    }

    os << "}";

    if (sequence_ids.empty()) {
        os << "}";
        os.close();
        return;
    }

    os << ",\"piles\":{";
    is_first = true;
    for (const auto& it: sequence_ids) {
        if (!is_first) {
            os << ",";
        }
        is_first = false;

        os << piles_[it]->to_json();
    }

    os << "}}";
    os.close();
}

void Graph::print_debug(const std::string& prefix) const {
    if (!prefix.empty()) {
        print_csv(prefix + ".csv");
        print_json(prefix + ".json");
    }
}

}
