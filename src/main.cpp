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

    std::vector<bool> is_valid_read, is_valid_overlap;
    prefilterData(is_valid_read, is_valid_overlap, reads_path, overlaps_path, overlap_type);

    std::vector<std::shared_ptr<Read>> reads;
    std::vector<std::shared_ptr<Overlap>> overlaps;
    preprocessData(reads, overlaps, is_valid_read, is_valid_overlap, reads_path,
        overlaps_path, overlap_type);

    auto graph = createGraph(reads, overlaps);
    overlaps.clear();
    reads.clear();

    graph->remove_isolated_nodes();
    graph->remove_transitive_edges();
    // graph->create_unitigs();
    // graph->remove_long_edges();
    uint32_t r = 0;
    fprintf(stderr, "\n");
    while (r < 10) {
        fprintf(stderr, "Simplification round %d {\n", r);
        graph->create_unitigs();
        graph->remove_tips();
        graph->remove_bubbles();
        ++r;
        fprintf(stderr, "}\n\n");
    }
    graph->print_contigs();

    //graph->print_dot();
    //graph->print_csv();

    return 0;
}
