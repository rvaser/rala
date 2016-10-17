#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include <unordered_set>
#include <algorithm>

#include "overlap.hpp"
#include "graph.hpp"
#include "read.hpp"
#include "bioparser/src/bioparser.hpp"

using namespace RALA;

int main(int argc, char** argv) {

    std::vector<std::shared_ptr<Overlap>> overlaps;
    auto oreader = BIOPARSER::createReader<Overlap, BIOPARSER::MhapReader>(argv[1]);
    oreader->read_objects(overlaps, 1000000000);

    std::vector<std::shared_ptr<Read>> reads;
    auto rreader = BIOPARSER::createReader<Read, BIOPARSER::FastqReader>(argv[2]);
    rreader->read_objects(reads, 1000000000);

    preprocessData(reads, overlaps);

    auto graph = createGraph(reads, overlaps);
    for (auto& read: reads) {
        if (read != nullptr) read.reset();
    }
    for (auto& overlap: overlaps) {
        if (overlap != nullptr) overlap.reset();
    }

    graph->remove_isolated_nodes();
    graph->remove_transitive_edges();
    graph->create_unitigs();
    //graph->remove_long_edges();
    uint32_t r = 0;
    while (r < 10) {
        //graph->create_unitigs();
        graph->remove_tips();
        graph->remove_bubbles();
        ++r;
    }
    //graph->create_unitigs();
    //graph->remove_tips();
    //graph->create_unitigs();
    //graph->print_contigs();

    //graph->print_dot();
    graph->print_csv();

    return 0;
}
