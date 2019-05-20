#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <algorithm>

#include "sequence.hpp"
#include "bioparser/bioparser.hpp"

static const char* version = "v0.7.0";

static struct option options[] = {
    {"reference", required_argument, 0, 'r'},
    {"reference-length", required_argument, 0, 'l'},
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    std::string reference_path = "";
    uint32_t reference_length = 0;

    char opt;
    while ((opt = getopt_long(argc, argv, "r:l:h", options, nullptr)) != -1) {
        switch (opt) {
            case 'r':
                reference_path = optarg;
                break;
            case 'l':
                reference_length = atoi(optarg);
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

    if (input_paths.size() == 0) {
        fprintf(stderr, "[rast::] error: missing input file(s)!\n");
        help();
        exit(1);
    }

    auto calculate_statistics = [](std::vector<std::unique_ptr<rala::Sequence>>& src,
        uint32_t length) -> uint32_t {

        std::sort(src.begin(), src.end(),
            [](const std::unique_ptr<rala::Sequence>& lhs,
                const std::unique_ptr<rala::Sequence>& rhs) {
                return lhs->data().size() > rhs->data().size();
            });

        uint32_t total_length = 0;
        for (const auto& it: src) {
            total_length += it->data().size();
            if (total_length > length / 2) {
                return it->data().size();
            }
        }

        return 0;
    };

    auto is_suffix = [](const std::string& src, const std::string& suffix) -> bool {
        if (src.size() < suffix.size()) {
            return false;
        }
        return src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    uint32_t reference_ng50 = 0;

    if (!reference_path.empty()) {

        std::unique_ptr<bioparser::Parser<rala::Sequence>> rparser = nullptr;

        if (is_suffix(reference_path, ".fasta") ||
            is_suffix(reference_path, ".fa") ||
            is_suffix(reference_path, ".fasta.gz") ||
            is_suffix(reference_path, ".fa.gz")) {

            rparser = bioparser::createParser<bioparser::FastaParser,
                rala::Sequence>(reference_path);

        } else if (is_suffix(reference_path, ".fastq") ||
            is_suffix(reference_path, ".fq") ||
            is_suffix(reference_path, ".fastq.gz") ||
            is_suffix(reference_path, ".fq.gz")) {

            rparser = bioparser::createParser<bioparser::FastqParser,
                rala::Sequence>(reference_path);

        } else {
            fprintf(stderr, "[rast::] error: file %s has unsupported format "
                "extension (valid extensions: .fasta, .fasta.gz, .fa, .fa.gz, "
                ".fastq, .fastq.gz, .fq, .fq.gz)!\n", reference_path.c_str());
            exit(1);
        }

        std::vector<std::unique_ptr<rala::Sequence>> reference;
        rparser->parse(reference, -1);

        reference_length = 0;
        for (const auto& it: reference) {
            reference_length += it->data().size();
        }
        reference_ng50 = calculate_statistics(reference, reference_length);

        fprintf(stderr, "[rast::] number of reference sequences = %zu\n", reference.size());
        fprintf(stderr, "[rast::] reference length = %u\n", reference_length);
        fprintf(stderr, "[rast::] reference NG50 = %u\n", reference_ng50);
    }

    for (const auto& it: input_paths) {

        std::unique_ptr<bioparser::Parser<rala::Sequence>> aparser = nullptr;

        if (is_suffix(it, ".fasta") || is_suffix(it, ".fa") ||
            is_suffix(it, ".fasta.gz") || is_suffix(it, ".fa.gz")) {

            aparser = bioparser::createParser<bioparser::FastaParser,
                rala::Sequence>(it);

        } else if (is_suffix(it, ".fastq") || is_suffix(it, ".fq") ||
            is_suffix(it, ".fastq.gz") || is_suffix(it, ".fq.gz")) {

                aparser = bioparser::createParser<bioparser::FastqParser,
                    rala::Sequence>(it);

        } else {
            fprintf(stderr, "[rast::] warning: file %s has unsupported format "
                "extension (valid extensions: .fasta, .fasta.gz, .fa, .fa.gz, "
                ".fastq, .fastq.gz, .fq, .fq.gz)!\n", it.c_str());
            continue;
        }

        std::string base_name = it.substr(it.rfind('/') + 1);

        fprintf(stderr, "[rast::] processing assembly %s\n", base_name.c_str());

        std::vector<std::unique_ptr<rala::Sequence>> assembly;
        aparser->parse(assembly, -1);

        uint32_t assembly_length = 0;
        for (const auto& it: assembly) {
            assembly_length += it->data().size();
        }

        uint32_t assembly_n50 = calculate_statistics(assembly, assembly_length);

        fprintf(stderr, "[rast::] number of assembly sequences = %zu\n", assembly.size());
        fprintf(stderr, "[rast::] assembly length = %u\n", assembly_length);
        fprintf(stderr, "[rast::] assembly N50 = %u\n", assembly_n50);

        if (reference_length != 0) {
            uint32_t assembly_ng50 = calculate_statistics(assembly, reference_length);
            fprintf(stderr, "[rast::] assembly NG50 = %u\n", assembly_ng50);
            if (reference_ng50 != 0) {
                fprintf(stderr, "[rast::] assembly status = [%s]\n",
                    assembly_ng50 > 0.9 * reference_ng50 ? "OK" : "FAILED");
            }
        }
    }

    return 0;
}

void help() {
    printf(
        "usage: rast [options ...] <assembly> [<assembly> ...]\n"
        "\n"
        "    <assembly>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing an assembly\n"
        "\n"
        "    options:\n"
        "        -r, --reference <string>\n"
        "            file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "            containing the reference genome\n"
        "        -l, --reference-length <int>\n"
        "            approximate length of the genome if the reference is unknown\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n");
}
