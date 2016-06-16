/*!
 * @file Graph.hpp
 *
 * @brief Graph class header file
 */

#pragma once

#include <memory>
#include <vector>

namespace RALAY {

class Graph;
class Overlap;
std::unique_ptr<Graph> createGraph(const std::vector<std::unique_ptr<Overlap>>& overlaps);

class Graph {
public:

    ~Graph();

    void print() const;

    friend std::unique_ptr<Graph> createGraph(const std::vector<std::unique_ptr<Overlap>>& overlaps);

private:

    Graph(const std::vector<std::unique_ptr<Overlap>>& overlaps);
    Graph(const Graph&) = delete;
    const Graph& operator=(const Graph&) = delete;
};

}
