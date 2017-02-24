#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include <unordered_set>
#include <algorithm>

#include "overlap.hpp"
#include "graph.hpp"
#include "read.hpp"
#include "utils.hpp"
#include "bioparser/src/bioparser.hpp"

using namespace rala;

int main(int argc, char** argv) {

    std::string reads_path = argc > 2 ? argv[2] : "";
    std::string overlaps_path = argc > 3 ? argv[3] : "";
    uint32_t overlap_type = argc > 4 ? atoi(argv[4]) : 0;

    switch (atoi(argv[1])) {
        case 1:
            findChimericReads(reads_path, overlaps_path, overlap_type);
            return 0;
        case 2:
            findUncontainedReads(reads_path, overlaps_path, overlap_type);
            return 0;
        case 3:
            fastaToFastq(reads_path);
            return 0;
        case 0:
        default:
            break;
    }

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

    graph->print_csv("layout_graph.csv");

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

    return 0;
}
