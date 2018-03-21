#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "sequence.hpp"
#include "graph.hpp"
#include "thread_pool/thread_pool.hpp"

static const char* version = "v0.7.0";

static struct option options[] = {
    {"threads", required_argument, 0, 't'},
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    uint32_t num_threads = 1;

    char opt;
    while ((opt = getopt_long(argc, argv, "t:h", options, nullptr)) != -1) {
        switch (opt) {
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'v':
                printf("%s\n", version);
                exit(0);
            case 'h':
                help();
                exit(0);
            default:
                exit(1);
        }
    }

    std::vector<std::string> input_paths;

    for (int32_t i = optind; i < argc; ++i) {
        input_paths.emplace_back(argv[i]);
    }

    if (input_paths.size() < 2) {
        fprintf(stderr, "[rala::] error: missing input file(s)!\n");
        help();
        exit(1);
    }

    auto graph = rala::createGraph(input_paths[0], input_paths[1], num_threads);
    graph->construct();
    graph->simplify();

    std::vector<std::unique_ptr<rala::Sequence>> contigs;
    graph->extract_contigs(contigs);

    graph->print_csv("assembly_graph.csv");
    graph->print_gfa("assembly_graph.gfa");

    for (const auto& it: contigs) {
        fprintf(stdout, "%s\n%s\n", it->name().c_str(), it->data().c_str());
    }

    return 0;
}

void help() {
    printf(
        "usage: rala [options ...] <sequences> <overlaps>\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing sequences\n"
        "    <overlaps>\n"
        "        input file in MHAP/PAF format (can be compressed with gzip)\n"
        "        containing pairwise overlaps\n"
        "\n"
        "    options:\n"
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n");
}
