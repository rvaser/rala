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
class Graph;
std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps);

class Node;
class Edge;
class Graph {
public:

    ~Graph();

    void print() const;

    friend std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
        const std::vector<std::shared_ptr<Overlap>>& overlaps);

private:

    Graph(const std::vector<std::shared_ptr<Read>>& reads,
        const std::vector<std::shared_ptr<Overlap>>& overlaps);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;

    std::vector<std::unique_ptr<Node>> nodes_;
    std::vector<std::unique_ptr<Edge>> edges_;
};

}
