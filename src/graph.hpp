/*!
 * @file Graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <memory>
#include <vector>
#include <unordered_set>

namespace RALA {

class Read;
class ReadFrame;
class Overlap;

/*
 * @brief Temporary function for chimeric read detection with reference
 */
void findChimeriReads(const std::string& overlaps_path, uint32_t overlap_type);

/*
 * @brief Remove contained reads and their overlaps; Remove unreliable overlaps
 */
void prefilterData(std::vector<bool>& is_valid_read, std::vector<bool>& is_valid_overlap,
    const std::string& reads_path, const std::string& overlaps_path, uint32_t overlap_type);

/*
 * @brief Trims reads from both sides based on read to read overlaps (Li 2016);
 */
void preprocessData(std::vector<std::shared_ptr<Read>>& reads, std::vector<std::shared_ptr<Overlap>>& overlaps,
    std::vector<bool>& is_valid_read, std::vector<bool>& is_valid_overlap, const std::string& reads_path,
    const std::string& overlaps_path, uint32_t overlap_type);

class Graph;
std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps);

class Graph {
public:

    ~Graph();

    /*
     * @brief Removes nodes without edges (no information loss)
     */
    void remove_isolated_nodes();

    /*
     * @brief Removes transitive edge (no information loss)
     * (inspired by Myers 1995 & 2005)
     */
    void remove_transitive_edges();

    /*
     * @brief Removes long edges (i.e. small overlaps, possible information loss and graph fragmentation)
     * (Li 2016)
     */
    void remove_long_edges();

    /*
     * @brief Removes tips (i.e. nodes with in_degree == 0 || out_degree == 0)
     */
    void remove_tips();

    /*
     * @brief Removes cycles (possible information loss)
     * (Tarjan 1972)
     */
    void remove_cycles();

    /*
     * @brief Removes bubbles (possible information loss and graph fragmentation)
     */
    uint32_t remove_bubbles();

    /*
     * @brief Creates unitigs by merging chains of overlapping reads
     */
    uint32_t create_unitigs();

    /*
     * @brief Outputs unitigs in FASTA format
     */
    void print_contigs() const;

    /*
     * @brief Prints assembly graph in graphviz format
     */
    void print_dot() const;

    /*
     * @brief Prints assembly graph in csv format
     */
    void print_csv() const;

    /*
     * @brief Temporary
     */
    void remove_selected_nodes_and_edges();

    friend std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
        const std::vector<std::shared_ptr<Overlap>>& overlaps);

private:

    Graph(const std::vector<std::shared_ptr<Read>>& reads,
        const std::vector<std::shared_ptr<Overlap>>& overlaps);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    void remove_marked_edges();

    void locate_bubble_sources(std::vector<uint32_t>& dst);

    class Node;
    class Edge;

    std::vector<std::shared_ptr<Node>> nodes_;
    std::vector<std::shared_ptr<Edge>> edges_;
    std::unordered_set<uint32_t> marked_edges_;
};

}
