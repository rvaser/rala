#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "graph.hpp"
#include "thread_pool/thread_pool.hpp"

static struct option options[] = {
    {"reads", required_argument, 0, 'i'},
    {"overlaps", required_argument, 0, 'j'},
    {"threads", required_argument, 0, 't'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    std::string reads_path;
    std::string overlaps_path;
    std::string output_path;
    uint32_t num_threads = std::thread::hardware_concurrency() / 2;

    while (true) {
        auto argument = getopt_long(argc, argv, "i:j:t:h", options, nullptr);
        if (argument == -1) {
            break;
        }

        switch (argument) {
            case 'i':
                reads_path = optarg;
                break;
            case 'j':
                overlaps_path = optarg;
                break;
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'h':
            default:
                help();
                return -1;
        }
    }

    if (reads_path.empty()) {
        fprintf(stderr, "rala:: error: missing option -i (reads file)!\n");
        help();
        return -1;
    }
    if (overlaps_path.empty()) {
        fprintf(stderr, "rala:: error: missing option -j (overlaps file)!\n");
        help();
        return -1;
    }

    auto graph = rala::createGraph(reads_path, overlaps_path, num_threads);
    graph->construct();

    return 0;
}

void help() {
    printf(
        "usage: rala -i <reads file> -j <overlaps file> [arguments ...]\n"
        "\n"
        "arguments:\n"
        "    -i, --reads <file>\n"
        "        (required)\n"
        "        input FASTQ file containg reads\n"
        "    -j, --overlaps <file>\n"
        "        (required)\n"
        "        input PAF file containing pairwise overlaps\n"
        "    -t, --threads <int>\n"
        "        default: hardware concurrency / 2\n"
        "        number of threads\n"
        "    -h, --help\n"
        "        prints out the help\n");
}
