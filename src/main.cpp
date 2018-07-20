#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "sequence.hpp"
#include "graph.hpp"
#include "thread_pool/thread_pool.hpp"

static const char* version = "v0.7.0";

static struct option options[] = {
    {"include-unassembled", no_argument, 0, 'u'},
    {"debug", required_argument, 0, 'd'},
    {"repeat-overlaps", required_argument, 0, 'r'},
    {"threads", required_argument, 0, 't'},
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    uint32_t num_threads = 1;
    bool drop_unassembled_sequences = true;
    std::string debug_prefix = "";
    std::string repeat_overlaps_path = "";

    char opt;
    while ((opt = getopt_long(argc, argv, "ud:r:t:h", options, nullptr)) != -1) {
        switch (opt) {
            case 'u':
                drop_unassembled_sequences = false;
                break;
            case 'd':
                debug_prefix = optarg;
                break;
            case 'r':
                repeat_overlaps_path = optarg;
                break;
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
    graph->construct(repeat_overlaps_path);
    graph->simplify();
    graph->print_debug(debug_prefix);

    std::vector<std::unique_ptr<rala::Sequence>> contigs;
    graph->extract_contigs(contigs, drop_unassembled_sequences);

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
        "        -u, --include-unassembled\n"
        "            output unassembled sequences (singletons and short contigs)\n"
        "        -d, --debug <string>\n"
        "            enable debug output with given prefix\n"
        "        -r, --repeat-overlaps <file>\n"
        "            input file in MHAP/PAF format (can be compress with gzip)\n"
        "            containing overlaps to repeat reads\n"
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n");
}
