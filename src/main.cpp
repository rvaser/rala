#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "read.hpp"
#include "graph.hpp"
#include "thread_pool/thread_pool.hpp"

static struct option options[] = {
    {"threads", required_argument, 0, 't'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    uint32_t num_threads = std::thread::hardware_concurrency() / 2;

    char opt;
    while ((opt = getopt_long(argc, argv, "t:h", options, nullptr)) != -1) {
        switch (opt) {
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'h':
            default:
                help();
                return -1;
        }
    }

    std::vector<std::string> input_paths;

    for (int32_t i = optind; i < argc; ++i) {
        input_paths.emplace_back(argv[i]);
    }

    if (input_paths.size() < 2) {
        fprintf(stderr, "rala:: error: missing input file(s)!\n");
        help();
        return -1;
    }

    auto graph = rala::createGraph(input_paths[0], input_paths[1], num_threads);
    graph->construct();
    graph->simplify();

    std::vector<std::unique_ptr<rala::Read>> contigs;
    graph->get_contigs(contigs);

    for (const auto& it: contigs) {
        fprintf(stdout, "%s\n%s\n", it->name().c_str(), it->sequence().c_str());
    }

    return 0;
}

void help() {
    printf(
        "usage: rala [options ...] <reads> <overlaps>\n"
        "\n"
        "    <reads>\n"
        "        input file in FASTA/FASTQ format containing reads\n"
        "    <overlaps>\n"
        "        input file in MHAP/PAF format containing pairwise overlaps\n"
        "        !note: if you are using an overlapper with the PAF file format,\n"
        "            reformat the read set with misc/fasta_formatter.py (or\n"
        "            misc/fastq_formatter.py) before running the overlapper\n"
        "\n"
        "    options:\n"
        "        -t, --threads <int>\n"
        "            default: hardware concurrency / 2\n"
        "            number of threads\n"
        "        -h, --help\n"
        "            prints out the help\n");
}
