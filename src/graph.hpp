/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <memory>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace bioparser {
    template<class T>
    class Parser;
}

namespace thread_pool {
    class ThreadPool;
}

namespace rala {

class Sequence;
class Pile;
class Overlap;

class Graph;
std::unique_ptr<Graph> createGraph(const std::string& sequences_path,
    const std::string& overlaps_path, uint32_t num_threads);

class Graph {
public:
    ~Graph();

    /*!
     * @brief Constructs the assembly graph by removing contained sequences and
     * transitive overlaps (removes chimeric and repetitive sequences before
     * construction if flag is set)
     */
    void construct(bool preprocess = true);

    /*!
     * @brief Removes transitive edges and tips, pops bubbles
     */
    void simplify();

    /*!
     * @brief Removes nodes without edges (no information loss)
     */
    uint32_t remove_isolated_nodes();

    /*!
     * @brief Removes transitive edge (no information loss)
     * (inspired by Myers 1995 & 2005)
     */
    uint32_t remove_transitive_edges();

    /*!
     * @brief Removes long edges (i.e. small overlaps, possible information
     * loss and graph fragmentation) (Li 2016)
     */
    uint32_t remove_long_edges();

    /*!
     * @brief Removes nodes which are dead ends in graph
     */
    uint32_t remove_tips();

    /*!
     * @brief Removes bubbles
     */
    uint32_t remove_bubbles();

    /*!
     * @brief Creates unitigs by merging chains of overlapping sequences
     */
    uint32_t create_unitigs();

    /*!
     * @brief Stores all contigs into dst
     */
    void extract_contigs(std::vector<std::unique_ptr<Sequence>>& dst) const;

    /*!
     * @brief Prints assembly graph in csv format
     */
    void print_csv(std::string path) const;

    /*!
     * @brief Prints assembly graph in GFA format
     */
    void print_gfa(std::string path) const;

    /*!
     * @brief Prints all unresolved graph junctions into graphs/ folder if
     * such exists
     */
    void print_knots() const; // TODO: reimplement

    friend std::unique_ptr<Graph> createGraph(const std::string& sequences_path,
        const std::string& overlaps_path, uint32_t num_threads);
private:
    Graph(std::unique_ptr<bioparser::Parser<Sequence>> sparser,
        std::unique_ptr<bioparser::Parser<Overlap>> oparser,
        uint32_t num_threads);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    /*!
     * @brief Initializes all structures and trims sequences
     */
    void initialize();

    /*!
     * @brief Splits chimeric sequences and removes overlaps between sequences
     * that do not bridge repetitive genomic regions
     */
    void preprocess();

    uint64_t find_edge(uint64_t src, uint64_t dst);

    /*!
     * @brief Finds edges in path which do not affect the connectivity of the
     * rest of the graph if removed
     */
    void find_removable_edges(std::vector<uint64_t>& dst,
        const std::vector<uint64_t>& path);

    void remove_marked_objects();

    class Node;
    class Edge;

    std::unique_ptr<bioparser::Parser<Sequence>> sparser_;
    std::unordered_map<std::string, uint64_t> name_to_id_;

    std::vector<std::unique_ptr<Pile>> piles_;
    uint32_t coverage_median_;

    std::unique_ptr<bioparser::Parser<Overlap>> oparser_;
    std::vector<bool> is_valid_overlap_;

    std::unique_ptr<thread_pool::ThreadPool> thread_pool_;

    std::vector<std::unique_ptr<Node>> nodes_;
    std::unordered_set<uint64_t> marked_nodes_;

    std::vector<std::unique_ptr<Edge>> edges_;
    std::unordered_set<uint64_t> marked_edges_;
};

}
