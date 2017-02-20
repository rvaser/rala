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

    std::string reads_path = argv[1];
    std::string overlaps_path = argv[2];
    uint32_t overlap_type = atoi(argv[3]);

    // findChimericReads(reads_path, overlaps_path, overlap_type);
    // return 0;

    std::vector<std::shared_ptr<Read>> reads;
    std::vector<std::shared_ptr<Overlap>> overlaps;
    preprocessData(reads, overlaps, reads_path, overlaps_path, overlap_type);

    auto graph = createGraph(reads, overlaps);
    overlaps.clear();
    reads.clear();

    graph->remove_isolated_nodes();
    graph->remove_transitive_edges();

    while (true) {
        uint32_t num_changes = graph->remove_tips();
        num_changes += graph->remove_chimeras();
        num_changes += graph->remove_bubbles();
        num_changes += graph->create_unitigs();
        if (num_changes == 0) {
            break;
        }
    }

    // graph->remove_selected_nodes_and_edges();
    graph->remove_long_edges();
    while (true) {
        uint32_t num_changes = graph->create_unitigs();
        num_changes += graph->remove_tips();
        if (num_changes == 0) {
            break;
        }
    }
    graph->print_contigs();
    // graph->print_csv("layout_graph.csv");

    return 0;
}
