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
#include "timer.hpp"

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"

namespace rala {

constexpr uint32_t kChunkSize = 1024 * 1024 * 1024; // ~ 1GB

// assembly graph
constexpr double kTransitiveEdgeEps = 0.12;
constexpr uint32_t kMaxBubbleLength = 5000000;
constexpr uint32_t kMinUnitigSize = 6;
constexpr double kMaxOverlapRatio = 0.9;

bool comparable(double a, double b, double eps) {
    return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
        (b >= a * (1 - eps) && b <= a * (1 + eps));
}

template<class T>
void shrinkVector(std::vector<std::unique_ptr<T>>& src, uint64_t begin) {

    uint64_t i = begin, j = begin;
    for (; i < src.size(); ++i) {
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
            src[i].swap(src[j]);
        }
    }
    if (i < src.size()) {
        src.resize(i);
    }
}

class Graph::Node {
public:
    // Encapsulating read
    Node(uint32_t _id, uint32_t _read_id, const std::string& read_name,
        const std::string& _sequence)
            : id(_id), read_id(_read_id), pair(), sequence(_sequence),
            prefix_edges(), suffix_edges(), read_ids(1, _read_id),
            unitig_size(1), mark(false) {

        auto is_unitig = read_name.find("Utg=");
        if (is_unitig != std::string::npos) {
            unitig_size = atoi(read_name.c_str() + is_unitig + 4);
        }
    }
    // Unitig
    Node(uint32_t _id, Node* begin_node, Node* end_node,
        std::unordered_set<uint32_t>& marked_edges);
    // Circular unitig
    Node(uint32_t _id, Node* begin_node,
        std::unordered_set<uint32_t>& marked_edges);
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
        return (out_degree() > 0 && in_degree() == 0 &&
            unitig_size < kMinUnitigSize);
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
    // Encapsulating overlap
    Edge(uint32_t _id, Node* _begin_node, Node* _end_node, uint32_t _length)
            : id(_id), pair(), begin_node(_begin_node), end_node(_end_node),
            length(_length), mark(false) {
    }
    Edge(const Edge&) = delete;
    const Edge& operator=(const Edge&) = delete;

    ~Edge() {}

    std::string label() const {
        return begin_node->sequence.substr(0, length);
    }

    uint32_t id;
    Edge* pair;
    Node* begin_node;
    Node* end_node;
    uint32_t length;
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

std::unique_ptr<Graph> createGraph(const std::string& reads_path,
    const std::string& overlaps_path, uint32_t num_threads) {

    return std::unique_ptr<Graph>(new Graph(reads_path, overlaps_path,
        num_threads));
}

Graph::Graph(const std::string& reads_path, const std::string& overlaps_path,
    uint32_t num_threads)
        : rreader_(nullptr), read_infos_(), coverage_median_(0),
        oreader_(nullptr), overlap_infos_(), thread_pool_(nullptr), nodes_(),
        edges_(), marked_edges_() {

    thread_pool_ = thread_pool::createThreadPool(num_threads);

    auto extension = reads_path.substr(reads_path.rfind('.'));
    if (extension == ".fasta" || extension == ".fa") {
        rreader_ = bioparser::createReader<Read, bioparser::FastaReader>(
            reads_path);
    } else if (extension == ".fastq" || extension == ".fq") {
        rreader_ = bioparser::createReader<Read, bioparser::FastqReader>(
            reads_path);
    } else {
        fprintf(stderr, "rala::Graph::Graph error: "
            "file %s has unsupported format extension (valid extensions: "
            ".fasta, .fa, .fastq, .fq)!\n",
            reads_path.c_str());
        exit(1);
    }

    extension = overlaps_path.substr(overlaps_path.rfind('.'));
    if (extension == ".paf") {
        oreader_ = bioparser::createReader<Overlap, bioparser::PafReader>(
            overlaps_path);
    } else if (extension == ".mhap") {
        oreader_ = bioparser::createReader<Overlap, bioparser::MhapReader>(
            overlaps_path);
    } else {
        fprintf(stderr, "rala::Graph::Graph error: "
            "file %s has unsupported format extension (valid extensions: "
            ".paf, .mhap)!\n",
            overlaps_path.c_str());
        exit(1);
    }
}

Graph::~Graph() {
}

void Graph::initialize() {

    fprintf(stderr, "Graph::initialize {\n");
    Timer timer; timer.start();

    std::vector<std::unique_ptr<ReadInfo>>().swap(read_infos_);

    // store overlaps
    std::vector<std::shared_ptr<Overlap>> current_overlaps;
    std::vector<std::vector<uint32_t>> shrunken_overlaps;
    uint32_t num_duplicate_overlaps = 0;

    auto remove_duplicate_overlaps = [&]() -> void {
        for (uint32_t i = 0; i < current_overlaps.size(); ++i) {
            if (!overlap_infos_[current_overlaps[i]->id()]) {
                continue;
            }
            for (uint32_t j = i + 1; j < current_overlaps.size(); ++j) {
                if (!overlap_infos_[current_overlaps[j]->id()] ||
                    current_overlaps[i]->b_id() != current_overlaps[j]->b_id()) {
                    continue;
                }
                ++num_duplicate_overlaps;

                if (current_overlaps[i]->length() >
                    current_overlaps[j]->length()) {
                    overlap_infos_[current_overlaps[j]->id()] = false;
                } else {
                    overlap_infos_[current_overlaps[i]->id()] = false;
                    break;
                }
            }
        }
        return;
    };

    auto shrink_valid_overlaps = [&]() -> void {
        for (const auto& it: current_overlaps) {
            if (overlap_infos_[it->id()]) {
                shrunken_overlaps[it->a_id()].push_back(
                    (it->a_begin() + 1) << 1 | 0);
                shrunken_overlaps[it->a_id()].push_back(
                    (it->a_end() - 1) << 1 | 1);
                shrunken_overlaps[it->b_id()].push_back(
                    (it->b_begin() + 1) << 1 | 0);
                shrunken_overlaps[it->b_id()].push_back(
                    (it->b_end() - 1) << 1 | 1);
            }
        }
        current_overlaps.clear();
    };

    std::vector<uint32_t> read_lengths;
    uint32_t num_self_overlaps = 0;

    std::vector<std::future<void>> thread_futures;

    oreader_->rewind();
    while (true) {
        std::vector<std::shared_ptr<Overlap>> overlaps;
        auto status = oreader_->read_objects(overlaps, kChunkSize);

        overlap_infos_.resize(overlap_infos_.size() + overlaps.size(), true);

        for (const auto& it: overlaps) {
            uint32_t max_read_id = std::max(it->a_id(), it->b_id());
            if (shrunken_overlaps.size() <= max_read_id) {
                shrunken_overlaps.resize(max_read_id + 1);
                read_lengths.resize(max_read_id + 1);
            }

            if (read_infos_.size() <= max_read_id) {
                read_infos_.resize(max_read_id + 1);
            }

            read_lengths[it->a_id()] = it->a_length();

            // self overlap check TODO: remove chimeric reads!
            if (it->a_id() == it->b_id()) {
                overlap_infos_[it->id()] = false;
                ++num_self_overlaps;
                continue;
            }

            read_lengths[it->b_id()] = it->b_length();

            if (current_overlaps.size() != 0 &&
                current_overlaps.front()->a_id() != it->a_id()) {

                remove_duplicate_overlaps();
                shrink_valid_overlaps();
            }
            current_overlaps.push_back(it);
        }

        overlaps.clear();

        if (!status) {
            remove_duplicate_overlaps();
            shrink_valid_overlaps();
        }

        // create missing and update all coverage graphs
        for (uint32_t i = 0; i < shrunken_overlaps.size(); ++i) {
            if (shrunken_overlaps[i].empty()) {
                continue;
            }

            if (read_infos_[i] == nullptr) {
                read_infos_[i] = createReadInfo(i, read_lengths[i]);
            }

            thread_futures.emplace_back(thread_pool_->submit_task(
                [&](uint32_t id) -> void {
                    read_infos_[id]->update_coverage_graph(shrunken_overlaps[id]);
                    std::vector<uint32_t>().swap(shrunken_overlaps[id]);
                }, i));
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }
        thread_futures.clear();

        if (!status) {
            std::vector<std::vector<uint32_t>>().swap(shrunken_overlaps);
            break;
        }
    }

    // trim reads
    for (const auto& it: read_infos_) {
        if (it == nullptr) {
            continue;
        }

        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint32_t id) -> void {
                if (!read_infos_[id]->find_valid_region()) {
                    read_infos_[id].reset();
                };
            }, it->id()));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    uint32_t num_prefiltered_reads = 0;
    for (const auto& it: read_infos_) {
        if (it == nullptr) {
            ++num_prefiltered_reads;
        }
    }

    fprintf(stderr, "  number of self overlaps = %u\n", num_self_overlaps);
    fprintf(stderr, "  number of duplicate overlaps = %u\n", num_duplicate_overlaps);
    fprintf(stderr, "  number of prefiltered reads = %u\n", num_prefiltered_reads);
    timer.stop(); timer.print("  - time =");
    fprintf(stderr, "}\n");
}

void Graph::preprocess() {

    fprintf(stderr, "Graph::preprocess {\n");
    Timer timer; timer.start();

    // find coverage median of the dataset
    std::vector<std::future<void>> thread_futures;
    for (const auto& it: read_infos_) {
        if (it == nullptr) {
            continue;
        }

        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint32_t id) -> void {
                read_infos_[id]->find_coverage_median();
            }, it->id()));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    std::vector<uint16_t> medians;
    for (const auto& it: read_infos_) {
        if (it == nullptr) {
            continue;
        }
        medians.emplace_back(it->coverage_median());
    }

    std::sort(medians.begin(), medians.end());
    coverage_median_ = medians.size() % 2 == 1 ?
        medians[medians.size() / 2] :
        (medians[medians.size() / 2 - 1] + medians[medians.size() / 2]) / 2;
    fprintf(stderr, "  - dataset coverage median = %u\n", coverage_median_);

    // find chimeric reads
    for (const auto& it: read_infos_) {
        if (it == nullptr) {
            continue;
        }

        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint32_t id) -> void {
                if (!read_infos_[id]->find_coverage_pits(coverage_median_)) {
                    read_infos_[id].reset();
                }
            }, it->id()));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    // filter low quality reads
    uint32_t num_low_quality_reads = 0;
    for (auto& it: read_infos_) {
        if (it == nullptr) {
            continue;
        }
        if (it->coverage_median() * 2 < coverage_median_) {
            it.reset();
            ++num_low_quality_reads;
        }
    }
    fprintf(stderr, "  number of low quality reads = %u\n", num_low_quality_reads);

    // correct coverage graphs
    std::vector<std::vector<uint32_t>> shrunken_overlaps(read_infos_.size());
    oreader_->rewind();
    while (true) {
        std::vector<std::unique_ptr<Overlap>> overlaps;
        auto status = oreader_->read_objects(overlaps, kChunkSize);

        for (auto& it: overlaps) {
            if (!overlap_infos_[it->id()] ||
                read_infos_[it->a_id()] == nullptr ||
                read_infos_[it->b_id()] == nullptr) {
                it.reset();
                continue;
            }
            if(!it->update(read_infos_[it->a_id()]->begin(),
                read_infos_[it->a_id()]->end(),
                read_infos_[it->b_id()]->begin(),
                read_infos_[it->b_id()]->end())) {
                it.reset();
                continue;
            }

            shrunken_overlaps[it->a_id()].emplace_back(it->a_begin());
            shrunken_overlaps[it->a_id()].emplace_back(it->a_end());
            shrunken_overlaps[it->a_id()].emplace_back(it->b_id());
            shrunken_overlaps[it->a_id()].emplace_back(it->b_begin());
            shrunken_overlaps[it->a_id()].emplace_back(it->b_end());
            shrunken_overlaps[it->a_id()].emplace_back(it->orientation());

            shrunken_overlaps[it->b_id()].emplace_back(it->b_begin());
            shrunken_overlaps[it->b_id()].emplace_back(it->b_end());
            shrunken_overlaps[it->b_id()].emplace_back(it->a_id());
            shrunken_overlaps[it->b_id()].emplace_back(it->a_begin());
            shrunken_overlaps[it->b_id()].emplace_back(it->a_end());
            shrunken_overlaps[it->b_id()].emplace_back(it->orientation());

            it.reset();
        }

        for (uint32_t i = 0; i < shrunken_overlaps.size(); ++i) {
            if (shrunken_overlaps[i].empty()) {
                continue;
            }

            thread_futures.emplace_back(thread_pool_->submit_task(
                [&](uint32_t id) -> void {
                    read_infos_[id]->correct_coverage_graph(
                        shrunken_overlaps[id], read_infos_);
                    std::vector<uint32_t>().swap(shrunken_overlaps[id]);
                }, i));
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }
        thread_futures.clear();

        if (!status) {
            std::vector<std::vector<uint32_t>>().swap(shrunken_overlaps);

            // update coverage medians
            for (const auto& it: read_infos_) {
                if (it == nullptr) {
                    continue;
                }

                thread_futures.emplace_back(thread_pool_->submit_task(
                    [&](uint32_t id) -> void {
                        read_infos_[id]->find_coverage_median();
                    }, it->id()));
            }
            for (const auto& it: thread_futures) {
                it.wait();
            }
            thread_futures.clear();

            std::vector<uint16_t> medians;
            for (const auto& it: read_infos_) {
                if (it == nullptr) {
                    continue;
                }
                medians.emplace_back(it->coverage_median());
            }

            std::sort(medians.begin(), medians.end());
            coverage_median_ = medians.size() % 2 == 1 ?
                medians[medians.size() / 2] :
                (medians[medians.size() / 2 - 1] + medians[medians.size() / 2]) / 2;
            fprintf(stderr, "  - updated dataset coverage median = %u\n",
                coverage_median_);

            break;
        }
    }

    // find repetitive regions
    for (const auto& it: read_infos_) {
        if (it == nullptr) {
            continue;
        }

        thread_futures.emplace_back(thread_pool_->submit_task(
            [&](uint32_t id) -> void {
                read_infos_[id]->find_coverage_hills(coverage_median_);
            }, it->id()));
    }
    for (const auto& it: thread_futures) {
        it.wait();
    }
    thread_futures.clear();

    timer.stop(); timer.print("  - time =");
    fprintf(stderr, "}\n");
}

void Graph::construct(bool preprocess) {

    nodes_.clear();
    edges_.clear();
    marked_edges_.clear();

    this->initialize();
    if (preprocess) {
        this->preprocess();
    }

    fprintf(stderr, "Graph::construct {\n");
    Timer timer; timer.start();

    // store overlaps
    std::vector<std::unique_ptr<Overlap>> overlaps;
    oreader_->rewind();
    while (true) {
        uint64_t current_overlap_id = overlaps.size();
        auto status = oreader_->read_objects(overlaps, kChunkSize);

        for (uint64_t i = current_overlap_id; i < overlaps.size(); ++i) {
            auto& it = overlaps[i];
            if (!overlap_infos_[it->id()] ||
                read_infos_[it->a_id()] == nullptr ||
                read_infos_[it->b_id()] == nullptr) {

                it.reset();
                continue;
            }
            if(!it->update(read_infos_[it->a_id()]->begin(),
                read_infos_[it->a_id()]->end(),
                read_infos_[it->b_id()]->begin(),
                read_infos_[it->b_id()]->end())) {

                it.reset();
                continue;
            }

            switch (it->type()) {
                case OverlapType::kX:
                    it.reset();
                    break;
                case OverlapType::kB:
                    read_infos_[it->a_id()].reset();
                    it.reset();
                    break;
                case OverlapType::kA:
                    read_infos_[it->b_id()].reset();
                    it.reset();
                    break;
                default:
                    break;
            }

            if (it == nullptr) {
                continue;
            }

            // check for false overlaps
            if (!read_infos_[it->a_id()]->is_valid_overlap(it->a_begin(), it->a_end()) ||
                !read_infos_[it->b_id()]->is_valid_overlap(it->b_begin(), it->b_end())) {
                it.reset();
            }
        }

        shrinkVector(overlaps, current_overlap_id);

        if (!status) {
            // check if all non valid overlaps are deleted
            for (auto& it: overlaps) {
                if (read_infos_[it->a_id()] == nullptr ||
                    read_infos_[it->b_id()] == nullptr) {

                    it.reset();
                }
            }
            shrinkVector(overlaps, 0);

            break;
        }
    }
    std::vector<bool>().swap(overlap_infos_);

    // store reads
    std::vector<std::unique_ptr<Read>> reads;
    rreader_->rewind();
    while (true) {
        uint64_t current_read_id = reads.size();
        auto status = rreader_->read_objects(reads, kChunkSize);

        for (uint64_t i = current_read_id; i < reads.size(); ++i) {
            auto& it = reads[i];
            if (it->id() >= read_infos_.size()) {
                it.reset();
                continue;
            }
            if (read_infos_[it->id()] == nullptr) {
                it.reset();
                continue;
            }
            it->update(read_infos_[it->id()]->begin(),
                read_infos_[it->id()]->end());
            // read_infos_[it->id()].reset();
        }

        shrinkVector(reads, current_read_id);

        if (!status) {
            break;
        }
    }

    // create assembly graph
    std::vector<int32_t> read_id_to_node_id(reads.back()->id() + 1, -1);
    uint32_t node_id = 0;
    for (const auto& it: reads) {
        read_id_to_node_id[it->id()] = node_id;

        Node* node = new Node(node_id++, it->id(), it->name(), it->sequence());
        Node* node_complement = new Node(node_id++, it->id(), it->name(),
            it->reverse_complement());

        node->pair = node_complement;
        node_complement->pair = node;

        nodes_.push_back(std::unique_ptr<Node>(node));
        nodes_.push_back(std::unique_ptr<Node>(node_complement));
    }

    uint32_t edge_id = 0;
    for (const auto& it: overlaps) {

        auto node_a = nodes_[read_id_to_node_id[it->a_id()]].get();
        auto node_b = nodes_[read_id_to_node_id[it->b_id()] +
            it->orientation()].get();

        uint32_t a_begin = it->a_begin();
        uint32_t a_end = it->a_end();
        uint32_t b_begin = it->orientation() == 0 ? it->b_begin() :
            it->b_length() - it->b_end();
        uint32_t b_end = it->orientation() == 0 ? it->b_end() :
            it->b_length() - it->b_begin();

        if (it->type() == OverlapType::kAB) {
            Edge* edge = new Edge(edge_id++, node_a, node_b, a_begin - b_begin);
            Edge* edge_complement = new Edge(edge_id++, node_b->pair,
                node_a->pair, it->b_length() - b_end - it->a_length() + a_end);

            edge->pair = edge_complement;
            edge_complement->pair = edge;

            edges_.push_back(std::unique_ptr<Edge>(edge));
            edges_.push_back(std::unique_ptr<Edge>(edge_complement));

            node_a->suffix_edges.push_back(edge);
            node_a->pair->prefix_edges.push_back(edge_complement);
            node_b->prefix_edges.push_back(edge);
            node_b->pair->suffix_edges.push_back(edge_complement);

        } else if (it->type() == OverlapType::kBA) {
            Edge* edge = new Edge(edge_id++, node_b, node_a, b_begin - a_begin);
            Edge* edge_complement = new Edge(edge_id++, node_a->pair,
                node_b->pair, it->a_length() - a_end - it->b_length() + b_end);

            edge->pair = edge_complement;
            edge_complement->pair = edge;

            edges_.push_back(std::unique_ptr<Edge>(edge));
            edges_.push_back(std::unique_ptr<Edge>(edge_complement));

            node_b->suffix_edges.push_back(edge);
            node_b->pair->prefix_edges.push_back(edge_complement);
            node_a->prefix_edges.push_back(edge);
            node_a->pair->suffix_edges.push_back(edge_complement);
        }
    }

    // log
    fprintf(stderr, "  number of graph nodes = %zu\n", nodes_.size());
    fprintf(stderr, "  number of graph edges = %zu\n", edges_.size());
    timer.stop(); timer.print("  - time =");
    fprintf(stderr, "}\n");
}

void Graph::simplify() {

    fprintf(stderr, "Graph::simplify {\n");
    Timer timer; timer.start();

    this->remove_isolated_nodes();
    uint32_t num_transitive_edges = this->remove_transitive_edges();

    uint32_t num_tips = 0;
    uint32_t num_bubbles = 0;

    while (true) {
        uint32_t num_changes = this->create_unitigs();

        uint32_t num_changes_part = this->remove_tips();
        num_tips += num_changes_part;
        num_changes += num_changes_part;

        num_changes_part = this->remove_bubbles();
        num_bubbles += num_changes_part;
        num_changes += num_changes_part;

        if (num_changes == 0) {
            break;
        }
    }

    uint32_t num_long_edges = this->remove_long_edges();

    while (true) {
        uint32_t num_changes = this->create_unitigs();

        uint32_t num_changes_part = this->remove_tips();
        num_tips += num_changes_part;
        num_changes += num_changes_part;

        if (num_changes == 0) {
            break;
        }
    }

    fprintf(stderr, "  number of transitive edges = %u\n", num_transitive_edges);
    fprintf(stderr, "  number of tips = %u\n", num_tips);
    fprintf(stderr, "  number of bubbles = %u\n", num_bubbles);
    fprintf(stderr, "  number of long edges = %u\n", num_long_edges);

    timer.stop(); timer.print("  - time =");
    fprintf(stderr, "}\n");
}

uint32_t Graph::remove_isolated_nodes() {

    uint32_t num_isolated_nodes = 0;

    for (auto& it: nodes_) {
        if (it == nullptr) {
            continue;
        }
        if ((it->in_degree() == 0 && it->out_degree() == 0 &&
            it->unitig_size < kMinUnitigSize) || it->mark == true) {
            it.reset();
            ++num_isolated_nodes;
        }
    }

    return num_isolated_nodes;
}

uint32_t Graph::remove_transitive_edges() {

    uint32_t num_transitive_edges = 0;
    std::vector<Edge*> candidate_edge(nodes_.size(), nullptr);

    for (const auto& node_x: nodes_) {
        if (node_x == nullptr) {
            continue;
        }

        for (const auto& edge: node_x->suffix_edges) {
            candidate_edge[edge->end_node->id] = edge;
        }

        for (const auto& edge_xy: node_x->suffix_edges) {
            const auto& node_y = nodes_[edge_xy->end_node->id];

            for (const auto& edge_yz: node_y->suffix_edges) {
                uint32_t z = edge_yz->end_node->id;

                if (candidate_edge[z] != nullptr &&
                    candidate_edge[z]->mark == false) {

                    if (comparable(edge_xy->length + edge_yz->length,
                        candidate_edge[z]->length, kTransitiveEdgeEps)) {

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

    return num_transitive_edges;
}

uint32_t Graph::remove_long_edges() {

    uint32_t num_long_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || node->suffix_edges.size() < 2){
            continue;
        }

        for (const auto& edge1: node->suffix_edges) {
            for (const auto& edge2: node->suffix_edges) {
                if (edge1->id == edge2->id ||
                    edge1->mark == true ||
                    edge2->mark == true) {

                    continue;
                }
                if (node->length() - edge2->length <
                    (node->length() - edge1->length) * kMaxOverlapRatio) {

                    edge2->mark = true;
                    edge2->pair->mark = true;
                    marked_edges_.insert(edge2->id);
                    marked_edges_.insert(edge2->pair->id);
                    ++num_long_edges;
                }
            }
        }
    }
    remove_marked_edges();

    return num_long_edges;
}

uint32_t Graph::remove_tips() {

    uint32_t num_tip_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || !node->is_tip()) {
            continue;
        }

        uint32_t num_removed_edges = 0;

        for (const auto& edge: node->suffix_edges) {
            if (edge->end_node->in_degree() > 1) {
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

    return num_tip_edges;
}

uint32_t Graph::remove_bubbles() {

    std::vector<uint32_t> distance(nodes_.size(), 0);
    std::vector<uint32_t> visited(nodes_.size(), 0);
    uint32_t visited_length = 0;
    std::vector<int32_t> predecessor(nodes_.size(), -1);
    std::deque<uint32_t> queue;

    auto extract_path = [&](std::vector<int32_t>& dst, int32_t source,
        int32_t sink) -> void {

        int32_t curr_id = sink;
        while (curr_id != source) {
            dst.push_back(curr_id);
            curr_id = predecessor[curr_id];
        }
        dst.push_back(source);
        std::reverse(dst.begin(), dst.end());
    };

    auto calculate_path_length = [&](const std::vector<int32_t>& path)
        -> uint32_t {

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

    auto is_valid_bubble = [&](const std::vector<int32_t>& path,
        const std::vector<int32_t>& other_path) -> bool {

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
        if (std::min(path_length, other_path_length) <
            std::max(path_length, other_path_length) * 0.8) {

            for (uint32_t i = 1; i < other_path.size() - 1; ++i) {
                if (nodes_[other_path[i]]->in_degree() > 1 ||
                    nodes_[other_path[i]]->out_degree() > 1) {
                    return false;
                }
            }
            for (uint32_t i = 1; i < path.size() - 1; ++i) {
                if (nodes_[path[i]]->in_degree() > 1 ||
                    nodes_[path[i]]->out_degree() > 1) {
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
            std::vector<int32_t> path, other_path;
            extract_path(path, source, sink);
            other_path.push_back(sink);
            extract_path(other_path, source, sink_other_predecesor);

            if (is_valid_bubble(path, other_path)) {
                uint32_t path_num_reads = 0;
                for (const auto& it: path) {
                    path_num_reads += nodes_[it]->unitig_size;
                }

                uint32_t other_path_num_reads = 0;
                for (const auto& it: other_path) {
                    other_path_num_reads += nodes_[it]->unitig_size;
                }

                std::vector<uint32_t> edges_for_removal;
                if (path_num_reads > other_path_num_reads) {
                    find_removable_edges(edges_for_removal, other_path);
                } else {
                    find_removable_edges(edges_for_removal, path);
                }

                for (const auto& edge_id: edges_for_removal) {
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

    return num_bubbles_popped;
}

uint32_t Graph::create_unitigs() {

    uint32_t node_id = nodes_.size();
    std::vector<bool> visited(nodes_.size(), false);
    std::vector<std::unique_ptr<Node>> new_nodes;

    uint32_t num_unitigs_created = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || visited[node->id] || node->is_junction()) {
            continue;
        }

        bool is_circular = false;
        auto first_node = node.get();
        while (!first_node->is_junction()) {
            visited[first_node->id] = true;
            visited[first_node->pair->id] = true;
            if (first_node->in_degree() == 0 ||
                first_node->prefix_edges.front()->begin_node->is_junction()) {
                break;
            }
            first_node = first_node->prefix_edges.front()->begin_node;
            if (first_node->id == node->id) {
                is_circular = true;
                break;
            }
        }

        auto last_node = node.get();
        while (!last_node->is_junction()) {
            visited[last_node->id] = true;
            visited[last_node->pair->id] = true;
            if (last_node->out_degree() == 0 ||
                last_node->suffix_edges.front()->end_node->is_junction()) {
                break;
            }
            last_node = last_node->suffix_edges.front()->end_node;
            if (last_node->id == node->id) {
                is_circular = true;
                break;
            }
        }

        Node* unitig = nullptr;
        Node* unitig_complement = nullptr;

        if (is_circular) {
            unitig = new Node(node_id++, first_node, marked_edges_);
            unitig_complement = new Node(node_id++, first_node->pair,
                marked_edges_);
        } else if (first_node->id != last_node->id) {
            unitig = new Node(node_id++, first_node, last_node, marked_edges_);
            unitig_complement = new Node(node_id++, last_node->pair,
                first_node->pair, marked_edges_);
        }

        if (unitig != nullptr && unitig_complement != nullptr) {
            unitig->pair = unitig_complement;
            unitig_complement->pair = unitig;

            new_nodes.push_back(std::unique_ptr<Node>(unitig));
            new_nodes.push_back(std::unique_ptr<Node>(unitig_complement));

            ++num_unitigs_created;
        }
    }

    for (auto& node: new_nodes) {
        nodes_.push_back(std::move(node));
    }

    remove_marked_edges();
    remove_isolated_nodes();

    return num_unitigs_created;
}

void Graph::print_contigs() const {

    fprintf(stderr, "Graph::print_contigs {\n");

    uint32_t contig_id = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr || node->id % 2 == 0 ||
            node->unitig_size < kMinUnitigSize || node->length() < 10000) {
            continue;
        }
        fprintf(stderr, "  - contig %d, num reads = %u, length = %zu (%d -> %d)\n",
            contig_id, node->unitig_size, node->sequence.size(),
            node->read_ids.front(), node->read_ids.back());
        fprintf(stdout, ">Contig_%u_(Utg=%u:Len=%lu)\n%s\n", contig_id++,
            node->unitig_size, node->sequence.size(), node->sequence.c_str());
    }

    fprintf(stderr, "}\n");
}

int32_t Graph::find_edge(uint32_t src, uint32_t dst) {
    for (const auto& edge: nodes_[src]->suffix_edges) {
        if (edge->end_node->id == dst) {
            return edge->id;
        }
    }
    return -1;
}

void Graph::find_removable_edges(std::vector<uint32_t>& dst,
    const std::vector<int32_t>& path) {

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
        // remove whole path
        for (uint32_t i = 0; i < path.size() - 1; ++i) {
            dst.push_back(find_edge(path[i], path[i + 1]));
        }
        return;
    }

    if (pref != -1 && nodes_[path[pref]]->out_degree() > 1) {
        return;
    }
    if (suff != -1 && nodes_[path[suff]]->in_degree() > 1) {
        return;
    }

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
}

void Graph::print_csv(std::string path) const {

    auto graph_file = fopen(path.c_str(), "w");

    for (const auto& it: nodes_) {
        if (it == nullptr || it->id % 2 == 0) {
            continue;
        }
        fprintf(graph_file, "%u L:%u R:%d U:%d,%u L:%u R:%d U:%d,0,-\n",
            it->id, it->length(), it->read_id, it->unitig_size, it->pair->id,
            it->pair->length(), it->pair->read_id, it->pair->unitig_size);
    }

    for (const auto& it: edges_) {
        if (it == nullptr) {
            continue;
        }
        fprintf(graph_file, "%u L:%u R:%d U:%d,%u L:%u R:%d U:%d,1,%d %d\n",
            it->begin_node->id, it->begin_node->length(),
            it->begin_node->read_id, it->begin_node->unitig_size,
            it->end_node->id, it->end_node->length(), it->end_node->read_id,
            it->end_node->unitig_size, it->id, it->length);
    }

    fclose(graph_file);
}

void Graph::print_knots() const {

    std::vector<bool> visited(edges_.size(), false);

    for (const auto& node: nodes_) {
        if (node == nullptr || node->suffix_edges.size() < 2) {
            continue;
        }
        if (read_infos_[node->read_id] == nullptr) {
            continue;
        }

        std::vector<uint16_t> graph1(read_infos_[node->read_id]->coverage_graph());
        uint32_t begin1 = read_infos_[node->read_id]->begin();
        uint32_t end1 = read_infos_[node->read_id]->end();
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
            if (read_infos_[id] == nullptr) {
                continue;
            }

            std::vector<uint16_t> graph2(read_infos_[id]->coverage_graph());
            uint32_t begin2 = read_infos_[id]->begin();
            uint32_t end2 = read_infos_[id]->end();
            if (edge->end_node->id % 2 != 0) {
                std::reverse(graph2.begin(), graph2.end());
                uint32_t tmp = begin2;
                begin2 = graph2.size() - end2;
                end2 = graph2.size() - tmp;
            }

            std::ofstream out("graphs/e" + std::to_string(edge->id));
            out << "x " << node->read_id << " " << id << " median diff" << std::endl;
            for (uint32_t i = 0; i < begin1 + edge->length + (end2 - begin2); ++i) {
                uint32_t g1 = (i < graph1.size() ? graph1[i] : 0);
                uint32_t g2 = (i < begin1 + edge->length ? 0 :
                    graph2[i - (edge->length + begin1) + begin2]);
                out << i << " " << g1 << " " << g2 << " " << coverage_median_
                    << " " << (int32_t) (g1 - g2) << std::endl;
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
