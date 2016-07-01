/*!
 * @file Graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <memory>
#include <vector>

namespace RALAY {

class Read;
class Overlap;
/*
 * @brief Trims reads from both sides based on read to read overlaps
 */
void trimReads(std::vector<std::shared_ptr<Read>>& reads, std::vector<std::shared_ptr<Overlap>>& overlaps);

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
     * (inspired by Li 2016)
     */
    void remove_long_edges();

    /*
     * @brief Removes bubbles (possible information loss and graph fragmentation)
     * (inspired by Li 2016)
     */
    void remove_bubbles();

    /*
     * @brief Prints graph in graphviz format
     */
    void print() const;

    friend std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
        const std::vector<std::shared_ptr<Overlap>>& overlaps);

private:

    Graph(const std::vector<std::shared_ptr<Read>>& reads,
        const std::vector<std::shared_ptr<Overlap>>& overlaps);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    void remove_marked_edges();

    class Node;
    class Edge;

    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<std::unique_ptr<Edge>> edges_;
};

}
