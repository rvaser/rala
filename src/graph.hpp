/*!
 * @file Graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <memory>
#include <vector>

namespace RALA {

class Read;
class Overlap;
/*
 * @brief Trims reads from both sides based on read to read overlaps (taken from Miniasm);
 * Removes chimeric reads (taken from Miniasm)
 * (Li 2016)
 */
void preprocessData(std::vector<std::shared_ptr<Read>>& reads,
    std::vector<std::shared_ptr<Overlap>>& overlaps);

void calculateReadCoverages(std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps);

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
    void remove_bubbles();

    void print_contigs() const;
    void create_unitigs();

    /*
     * @brief Prints graph in graphviz format
     */
    void print_dot() const;

    /*
     * @brief Prints graph in csv format
     */
    void print_csv() const;

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

    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<std::unique_ptr<Edge>> edges_;
};

}
