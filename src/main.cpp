#include <stdio.h>
#include <stdlib.h>

#include "overlap.hpp"
#include "graph.hpp"
#include "read.hpp"
#include "preprocess.hpp"
#include "thread_pool/src/thread_pool.hpp"

using namespace rala;

int main(int argc, char** argv) {

    std::string reads_path = argc > 2 ? argv[2] : "";
    std::string overlaps_path = argc > 3 ? argv[3] : "";
    uint32_t overlap_type = argc > 4 ? atoi(argv[4]) : 0;
    std::string mappings_path = argc > 5 ? argv[5] : "";

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

    std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool();

    std::vector<std::shared_ptr<Read>> reads;
    std::vector<std::shared_ptr<Overlap>> overlaps;
    std::vector<std::vector<uint32_t>> mappings;
    double median;
    preprocessData(reads, overlaps, mappings, median, reads_path, overlaps_path,
        overlap_type, thread_pool);

    auto graph = createGraph(reads, overlaps);
    overlaps.clear();
    reads.clear();

    graph->remove_isolated_nodes();
    graph->remove_transitive_edges();

    graph->print_csv("layout_graph.csv");

    while (true) {
        uint32_t num_changes = graph->remove_tips();
        // num_changes += graph->remove_chimeras();
        num_changes += graph->remove_bubbles();
        // num_changes += graph->create_unitigs();
        if (num_changes == 0) {
            break;
        }
    }

    // graph->remove_selected_nodes_and_edges();
    // graph->print_knots(mappings, median);

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
