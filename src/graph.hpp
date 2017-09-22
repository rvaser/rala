/*!
 * @file graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>
#include <unordered_set>

namespace bioparser {
    template<class T>
    class Reader;
}

namespace thread_pool {
    class ThreadPool;
}

namespace rala {

class Read;
class ReadInfo;
class Overlap;

class Graph;
std::unique_ptr<Graph> createGraph(const std::string& reads_path,
    const std::string& overlaps_path, uint32_t num_threads);

class Graph {
public:
    ~Graph();

    /*!
     * @brief Constructs the assembly graph by removing contained reads and
     * transitive overlaps
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
     * @brief Creates unitigs by merging chains of overlapping reads
     */
    uint32_t create_unitigs();

    /*!
     * @brief Outputs unitigs in FASTA format
     */
    void print_contigs() const;

    void print_knots() const;

    /*!
     * @brief Prints assembly graph in csv format
     */
    void print_csv(std::string path) const;

    /*!
     * @brief For testing purposes
     */
    void remove_selected_nodes_and_edges();

    friend std::unique_ptr<Graph> createGraph(const std::string& reads_path,
        const std::string& overlaps_path, uint32_t num_threads);
private:
    Graph(const std::string& reads_path, const std::string& overlaps_path,
        uint32_t num_threads);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    /*!
     * @brief Initializes all structures by reading the overlaps and trimming
     * reads
     */
    void initialize();

    /*!
     * @brief Removes chimeric reads and those that do not bridge repetitive
     * genomic regions
     */
    void preprocess();

    int32_t find_edge(uint32_t src, uint32_t dst);

    /*!
     * @brief Finds edges in path which if removed do not affect the
     * connectivity of the rest of the graph
     */
    void find_removable_edges(std::vector<uint32_t>& dst,
        const std::vector<int32_t>& path);

    void remove_marked_edges();

    class Node;
    class Edge;

    std::unique_ptr<bioparser::Reader<Read>> rreader_;
    std::vector<std::unique_ptr<ReadInfo>> read_infos_;
    uint32_t coverage_median_;

    std::unique_ptr<bioparser::Reader<Overlap>> oreader_;
    std::vector<bool> overlap_infos_;

    std::unique_ptr<thread_pool::ThreadPool> thread_pool_;

    std::vector<std::shared_ptr<Node>> nodes_;
    std::vector<std::shared_ptr<Edge>> edges_;
    std::unordered_set<uint32_t> marked_edges_;
};

}
