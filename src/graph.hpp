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
     * @brief Simplifies graph without information loss
     * @details Removes nodes with in_degree = out_degree = 0 and transitive edges
     */
    void simplify();

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

    class Node;
    class Edge;

    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<std::unique_ptr<Edge>> edges_;
};

}
