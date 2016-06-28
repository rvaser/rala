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
    graph->print();

    return 0;
}
