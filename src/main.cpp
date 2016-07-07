#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include <unordered_set>

#include "overlap.hpp"
#include "graph.hpp"
#include "read.hpp"
#include "bioparser/src/bioparser.hpp"

using namespace RALAY;

int main(int argc, char** argv) {

    std::vector<std::shared_ptr<Overlap>> overlaps;
    auto reader = BIOPARSER::createMhapReader<Overlap>(argv[1]);
    reader->read_objects(overlaps, 1000000000);

    std::vector<std::shared_ptr<Read>> reads;
    auto qreader = BIOPARSER::createFastqReader<Read>(argv[2]);
    qreader->read_objects(reads, 1000000000);

    auto graph = createGraph(reads, overlaps);
    for (auto& it: reads) {
        it.reset();
    }
    for (auto& it: overlaps) {
        it.reset();
    }

    graph->remove_isolated_nodes();
    graph->remove_transitive_edges();
    //graph->remove_cycles();
    graph->remove_bubbles();
    fprintf(stderr, "ROUND 2\n");
    graph->remove_bubbles();
    fprintf(stderr, "ROUND 3\n");
    graph->remove_bubbles();
    fprintf(stderr, "ROUND 4\n");
    graph->remove_bubbles();
    fprintf(stderr, "ROUND 5\n");
    graph->remove_bubbles();
    graph->remove_tips();
    //graph->remove_long_edges();
    graph->create_unitigs();
    graph->remove_tips();
    graph->create_unitigs();
    graph->remove_tips();
    graph->create_unitigs();
    //graph->print_contigs();
    graph->print();

    return 0;
}
