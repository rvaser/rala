/*!
 * @file preprocess.cpp
 *
 * @brief Preprocessing source file
 */

#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <set>

#include "read.hpp"
#include "overlap.hpp"
#include "preprocess.hpp"
#include "bioparser/src/bioparser.hpp"
#include "thread_pool/src/thread_pool.hpp"

namespace rala {

constexpr uint32_t kChunkSize = 1024 * 1024 * 1024; // ~ 1GB

// preprocess
constexpr uint32_t kMinCoverage = 3;
constexpr double kMaxOverhangToOverlapRatio = 0.875;
constexpr double kSlopeRatio = 1.3;
constexpr uint32_t kSlopeWidth = 500;
constexpr double kSlopeWidthRatio = 0.05;
constexpr double kHillWidthRatio = 0.85;
constexpr double kMedianRatio = 2;

// unmapped reads
constexpr double kMinMappingPerc = 0.0;

template<class T>
void shrinkVector(std::vector<std::shared_ptr<T>>& dst, uint64_t start, const std::vector<bool>& is_valid) {

    uint64_t i = start, j = start;
    for (; i < dst.size(); ++i) {
        if (is_valid[dst[i]->id()]) {
            continue;
        }

        j = std::max(j, i);
        while (j < dst.size() && !is_valid[dst[j]->id()]) {
            ++j;
        }

        if (j >= dst.size()) {
            break;
        } else if (i != j) {
            dst[i].swap(dst[j]);
        }
    }
    if (i < dst.size()) {
        dst.resize(i);
    }
}

uint32_t removeDuplicateOverlaps(std::vector<bool>& is_valid_overlap,
    const std::vector<std::shared_ptr<Overlap>>& overlaps) {

    uint32_t num_duplicate_overlaps = 0;
    for (uint32_t i = 0; i < overlaps.size(); ++i) {
        if (is_valid_overlap[overlaps[i]->id()] == false) {
            continue;
        }
        for (uint32_t j = i + 1; j < overlaps.size(); ++j) {
            if (is_valid_overlap[overlaps[j]->id()] == false || overlaps[i]->b_id() != overlaps[j]->b_id()) {
                continue;
            }
            ++num_duplicate_overlaps;
            if (overlaps[i]->length() > overlaps[j]->length()) {
                is_valid_overlap[overlaps[j]->id()] = false;
            } else {
                is_valid_overlap[overlaps[i]->id()] = false;
                break;
            }
        }
    }

    return num_duplicate_overlaps;
};

uint32_t classifyOverlap(const std::shared_ptr<Overlap>& overlap) {

    uint32_t left_overhang = std::min(overlap->a_begin(), overlap->b_begin());
    uint32_t right_overhang = std::min(overlap->a_length() - overlap->a_end(),
        overlap->b_length() - overlap->b_end());

    uint32_t a_len = overlap->a_end() - overlap->a_begin();
    uint32_t b_len = overlap->b_end() - overlap->b_begin();

    if (a_len < (a_len + left_overhang + right_overhang) * kMaxOverhangToOverlapRatio ||
        b_len < (b_len + left_overhang + right_overhang) * kMaxOverhangToOverlapRatio) {
        return 0; // internal match
    }
    if (overlap->a_begin() <= overlap->b_begin() && (overlap->a_length() - overlap->a_end()) <= (overlap->b_length() - overlap->b_end())) {
        return 1; // a contained
    }
    if (overlap->a_begin() >= overlap->b_begin() && (overlap->a_length() - overlap->a_end()) >= (overlap->b_length() - overlap->b_end())) {
        return 2; // b contained
    }
    if (overlap->a_begin() > overlap->b_begin()) {
        return 3; // a to b overlap
    }
    return 4; // b to a overlap
}

void preprocessData(std::vector<std::shared_ptr<Read>>& reads, std::vector<std::shared_ptr<ReadInfo>>& read_infos,
    std::vector<std::shared_ptr<Overlap>>& overlaps, double& dataset_coverage_median,
    const std::string& reads_path, const std::string& overlaps_path, uint32_t overlap_type,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {

    fprintf(stderr, "Preprocessing data {\n");

    std::vector<std::shared_ptr<Overlap>> current_overlaps;
    std::vector<std::vector<uint32_t>> mappings;
    std::vector<bool> is_valid_overlap;
    uint32_t num_duplicate_overlaps = 0;
    uint32_t num_self_overlaps = 0;

    std::vector<uint32_t> read_lengths;

    auto oreader = overlap_type == 0 ?
        bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
        bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

    std::vector<std::shared_ptr<Overlap>> overlaps_of_spec;
    uint32_t spec = 70133;

    while (true) {
        auto status = oreader->read_objects(overlaps, kChunkSize);
        is_valid_overlap.resize(is_valid_overlap.size() + overlaps.size(), true);

        for (const auto& it: overlaps) {

            uint32_t max_read_id = std::max(it->a_id(), it->b_id());
            if (mappings.size() <= max_read_id) {
                mappings.resize(max_read_id + 1);
                read_lengths.resize(max_read_id + 1);
                read_infos.resize(max_read_id + 1);
            }

            read_lengths[it->a_id()] = it->a_length();
            read_lengths[it->b_id()] = it->b_length();

            // self overlap check
            if (it->a_id() == it->b_id()) {
                is_valid_overlap[it->id()] = false;
                ++num_self_overlaps;
                continue;
            }

            // duplicate overlaps check
            if (current_overlaps.size() != 0 && current_overlaps.front()->a_id() != it->a_id()) {
                num_duplicate_overlaps += removeDuplicateOverlaps(is_valid_overlap, current_overlaps);

                // save valid overlaps
                for (uint32_t i = 0; i < current_overlaps.size(); ++i) {
                    if (is_valid_overlap[current_overlaps[i]->id()]) {
                        const auto& valid_overlap = current_overlaps[i];

                        if (valid_overlap->a_id() == spec || valid_overlap->b_id() == spec) overlaps_of_spec.push_back(valid_overlap);

                        mappings[valid_overlap->a_id()].push_back(((valid_overlap->a_rc() ? valid_overlap->a_length() - valid_overlap->a_end() : valid_overlap->a_begin()) + 1) << 1 | 0);
                        mappings[valid_overlap->a_id()].push_back(((valid_overlap->a_rc() ? valid_overlap->a_length() - valid_overlap->a_begin() : valid_overlap->a_end()) - 1) << 1 | 1);

                        mappings[valid_overlap->b_id()].push_back(((valid_overlap->b_rc() ? valid_overlap->b_length() - valid_overlap->b_end() : valid_overlap->b_begin()) + 1) << 1 | 0);
                        mappings[valid_overlap->b_id()].push_back(((valid_overlap->b_rc() ? valid_overlap->b_length() - valid_overlap->b_begin() : valid_overlap->b_end()) - 1) << 1 | 1);
                    }
                }
                current_overlaps.clear();
            }
            current_overlaps.push_back(it);
        }

        overlaps.clear();

        // create missing coverage graphs
        {
            std::vector<std::future<std::unique_ptr<ReadInfo>>> thread_futures;
            std::vector<uint32_t> task_ids;
            for (uint32_t i = 0; i < mappings.size(); ++i) {
                if (read_infos[i] != nullptr) {
                    continue;
                }
                thread_futures.emplace_back(thread_pool->submit_task(
                    createReadInfo, i, read_lengths[i], std::ref(mappings[i])));
                task_ids.emplace_back(i);
            }

            for (uint32_t i = 0; i < thread_futures.size(); ++i) {
                thread_futures[i].wait();
                read_infos[task_ids[i]] = std::move(thread_futures[i].get());
                std::vector<uint32_t>().swap(mappings[task_ids[i]]);
            }
        }

        // update coverage graphs
        {
            std::vector<std::future<void>> thread_futures;
            std::vector<uint32_t> task_ids;
            for (uint32_t i = 0; i < mappings.size(); ++i) {
                if (read_infos[i] == nullptr) {
                    continue;
                }

                thread_futures.emplace_back(thread_pool->submit_task(
                    [&](uint32_t id) { read_infos[id]->update_coverage_graph(mappings[id]); }, i));
                task_ids.emplace_back(i);
            }

            for (uint32_t i = 0; i < thread_futures.size(); ++i) {
                thread_futures[i].wait();
                std::vector<uint32_t>().swap(mappings[task_ids[i]]);
            }
        }

        if (status == false) {
            break;
        }
    }
    oreader.reset();
    std::vector<std::vector<uint32_t>>().swap(mappings);

    // print overlaps of spec together
    /*{
        std::vector<std::vector<uint16_t>> coverage_graphs;
        coverage_graphs.push_back(read_infos[spec]->coverage_graph());
        for (const auto& it: overlaps_of_spec) {
            if (read_infos[it->a_id()] == nullptr || read_infos[it->b_id()] == nullptr) continue;

            uint32_t a_begin = it->a_rc() ? it->a_length() - it->a_end() : it->a_begin();
            uint32_t a_end = it->a_rc() ? it->a_length() - it->a_begin() : it->a_end();
            uint32_t b_begin = it->b_rc() ? it->b_length() - it->b_end() : it->b_begin();
            uint32_t b_end = it->b_rc() ? it->b_length() - it->b_begin() : it->b_end();

            uint32_t len = std::min(b_end - b_begin, a_end - a_begin);
            if (abs((a_end - a_begin) - (b_end - b_begin)) / (double) len > 0.01) continue;

            coverage_graphs.push_back(std::vector<uint16_t>(coverage_graphs.front().size(), 0));
            bool rc = it->a_rc() || it->b_rc();

            if (it->a_id() == spec) {
                std::vector<uint16_t> cg(read_infos[it->b_id()]->coverage_graph());
                for (uint32_t i = 0; i < len; ++i) {
                    if (rc == false) {
                        coverage_graphs.back()[a_begin + i] = cg[b_begin + i];
                    } else {
                        coverage_graphs.back()[a_begin + i] = cg[b_begin + len - 1 - i];
                    }
                }

            } else {
                std::vector<uint16_t> cg(read_infos[it->a_id()]->coverage_graph());
                for (uint32_t i = 0; i < len; ++i) {
                    if (rc == false) {
                        coverage_graphs.back()[b_begin + i] = cg[a_begin + i];
                    } else {
                        coverage_graphs.back()[b_begin + i] = cg[a_begin + len - 1 - i];
                    }
                }
            }
        }
        fprintf(stderr, "%d\n", coverage_graphs.size());
        std::ofstream out("graphs/h" + std::to_string(spec));
        for (uint32_t i = 0; i < coverage_graphs.front().size(); ++i) {
            out << i;
            for (uint32_t j = 0; j < coverage_graphs.size(); ++j) {
                out << " " << coverage_graphs[j][i];
            }
            out << std::endl;
        }
        out.close();
    }
    exit(1);*/

    fprintf(stderr, "  number of self overlaps = %u\n", num_self_overlaps);
    fprintf(stderr, "  number of duplicate overlaps = %u\n", num_duplicate_overlaps);

    // find longest contiguous read region which has coverage larger than kMinCoverage
    std::vector<bool> is_valid_read(read_infos.size(), true);
    {
        fprintf(stderr, "  Prefiltering data {\n");

        uint32_t num_prefiltered_reads = 0;
        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                is_valid_read[i] = false;
                ++num_prefiltered_reads;
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) { read_infos[id]->find_valid_region(kMinCoverage); }, i));
        }
        for (auto& it: thread_futures) {
            it.wait();
        }

        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            if (read_infos[i]->coverage_graph().empty()) {
                read_infos[i].reset();
                is_valid_read[i] = false;
                ++num_prefiltered_reads;
            }
        }

        fprintf(stderr, "    number of prefiltered reads = %u\n", num_prefiltered_reads);
        fprintf(stderr, "  }\n");
    }

    // find chimeric reads by looking for tight coverage pits
    {
        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) { read_infos[id]->find_coverage_pits(1.817, kSlopeWidth, kSlopeWidthRatio); }, i));
        }
        for (auto& it: thread_futures) {
            it.wait();
        }

        uint32_t num_chimeric_reads = 0;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }
            if (read_infos[i]->coverage_graph().empty()) {
                read_infos[i].reset();
                is_valid_read[i] = false;
                ++num_chimeric_reads;
            }
        }
        fprintf(stderr, "  number of chimeric reads = %u\n", num_chimeric_reads);
    }

    // correct coverage graphs
    {
        std::vector<std::shared_ptr<ReadInfo>> copies;
        for (const auto& it: read_infos) {
            copies.emplace_back(std::move(copyReadInfo(it)));
        }

        oreader = overlap_type == 0 ?
            bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
            bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

        while (true) {
            auto status = oreader->read_objects(overlaps, kChunkSize);

            for (const auto& it: overlaps) {
                if (!is_valid_overlap[it->id()] || !is_valid_read[it->a_id()] || !is_valid_read[it->b_id()]) {
                    continue;
                }

                bool is_valid = it->update(
                    read_infos[it->a_id()]->begin(),
                    read_infos[it->a_id()]->end(),
                    read_infos[it->b_id()]->begin(),
                    read_infos[it->b_id()]->end()
                );

                if (is_valid) {
                    uint32_t a_begin = it->a_begin() + (it->a_rc() ? read_infos[it->a_id()]->coverage_graph().size() - 1 - read_infos[it->a_id()]->end() : read_infos[it->a_id()]->begin());
                    uint32_t a_end = it->a_end() + (it->a_rc() ? read_infos[it->a_id()]->coverage_graph().size() - 1 -read_infos[it->a_id()]->end() : read_infos[it->a_id()]->begin());
                    if (it->a_rc()) {
                        auto tmp = a_begin;
                        a_begin = read_infos[it->a_id()]->coverage_graph().size() - 1 - a_end;
                        a_end = read_infos[it->a_id()]->coverage_graph().size() - 1 - tmp;
                    }
                    uint32_t b_begin = it->b_begin() + (it->b_rc() ? read_infos[it->b_id()]->coverage_graph().size() - 1 -read_infos[it->b_id()]->end() : read_infos[it->b_id()]->begin());
                    uint32_t b_end = it->b_end() + (it->b_rc() ? read_infos[it->b_id()]->coverage_graph().size() - 1 - read_infos[it->b_id()]->end() : read_infos[it->b_id()]->begin());
                    if (it->b_rc()) {
                        auto tmp = b_begin;
                        b_begin = read_infos[it->b_id()]->coverage_graph().size() - 1 - b_end;
                        b_end = read_infos[it->b_id()]->coverage_graph().size() - 1 - tmp;
                    }

                    copies[it->a_id()]->correct_coverage_graph(a_begin, a_end, read_infos[it->b_id()],
                        b_begin, b_end, it->a_rc() || it->b_rc());
                    copies[it->b_id()]->correct_coverage_graph(b_begin, b_end, read_infos[it->a_id()],
                        a_begin, a_end, it->a_rc() || it->b_rc());
                }
            }

            overlaps.clear();

            if (status == false) {
                break;
            }
        }

        for (uint32_t i = 0; i < copies.size(); ++i) {
            copies[i].swap(read_infos[i]);
        }

        oreader.reset();
    }

    // find coverage median and filter reads according to it
    {
        uint16_t median_of_medians = 0, average_median = 0;
        uint32_t num_low_median_reads = 0;

        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) { read_infos[id]->find_coverage_median(); }, i));
        }

        for (const auto& it: thread_futures) {
            it.wait();
        }

        std::vector<uint16_t> medians;
        uint32_t median_sum = 0;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            medians.emplace_back(read_infos[i]->coverage_median());
            median_sum += medians.back();
        }

        std::sort(medians.begin(), medians.end());
        median_of_medians = medians.size() % 2 == 1 ? medians[medians.size() / 2] :
            (medians[medians.size() / 2 - 1] + medians[medians.size() / 2]) / 2;
        average_median = median_sum / medians.size();

        fprintf(stderr, "  (median/average of coverage medians = %u/%u)\n", median_of_medians, average_median);

        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            if (read_infos[i]->coverage_median() * kMedianRatio < median_of_medians) {
                read_infos[i].reset();
                is_valid_read[i] = false;
                ++num_low_median_reads;
            }
        }

        fprintf(stderr, "  number of reads with coverage median < %u / %.3lf = %u\n",
            median_of_medians, kMedianRatio, num_low_median_reads);

        dataset_coverage_median = median_of_medians;
    }

    // smooth coverage graphs
    {
        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) { read_infos[id]->smooth_coverage_graph(); }, i));
        }
        for (auto& it: thread_futures) {
            it.wait();
        }
    }

    // find coverage hills which represent repetitive regions
    {
        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) { read_infos[id]->find_coverage_hills(kSlopeRatio, kSlopeWidth, kSlopeWidthRatio, kHillWidthRatio, dataset_coverage_median); }, i));
        }
        for (auto& it: thread_futures) {
            it.wait();
        }
    }

    std::vector<uint32_t> specials = {};
    for (const auto& id: specials) {
        if (read_infos[id] == nullptr) {
            continue;
        }
        read_infos[id]->print_csv("graphs/h" + std::to_string(id), dataset_coverage_median);
    }

    // reading valid overlaps into memory
    oreader = overlap_type == 0 ?
        bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
        bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

    while (true) {
        uint64_t current_overlap_id = overlaps.size();
        auto status = oreader->read_objects(overlaps, kChunkSize);

        for (uint64_t i = current_overlap_id; i < overlaps.size(); ++i) {
            auto& it = overlaps[i];
            if (!is_valid_overlap[it->id()]) {
                continue;
            }
            if (!is_valid_read[it->a_id()] || !is_valid_read[it->b_id()]) {
                is_valid_overlap[it->id()] = false;
                continue;
            }

            if (read_infos[it->a_id()]->coverage_hills().size() != 0) {
                uint32_t begin = it->a_rc() ? it->a_length() - it->a_end() : it->a_begin();
                uint32_t end = it->a_rc() ? it->a_length() - it->a_begin() : it->a_end();

                uint32_t valid_read_length = read_infos[it->a_id()]->end() - read_infos[it->a_id()]->begin();

                for (const auto& h: read_infos[it->a_id()]->coverage_hills()) {
                    if (begin < h.second && h.first < end) {
                        if (h.first < 0.10 * valid_read_length + read_infos[it->a_id()]->begin()) {
                            if (end < h.second) {
                                is_valid_overlap[it->id()] = false;
                            }
                        } else if (h.second > 0.9 * valid_read_length + read_infos[it->a_id()]->begin()) {
                            if (begin > h.first) {
                                is_valid_overlap[it->id()] = false;
                            }
                        }
                    }
                }
            }

            if (read_infos[it->b_id()]->coverage_hills().size() != 0) {
                uint32_t begin = it->b_rc() ? it->b_length() - it->b_end() : it->b_begin();
                uint32_t end = it->b_rc() ? it->b_length() - it->b_begin() : it->b_end();

                uint32_t valid_read_length = read_infos[it->b_id()]->end() - read_infos[it->b_id()]->begin();

                for (const auto& h: read_infos[it->b_id()]->coverage_hills()) {
                    if (begin < h.second && h.first < end) {
                        if (h.first < 0.10 * valid_read_length + read_infos[it->b_id()]->begin()) {
                            if (end < h.second) {
                                is_valid_overlap[it->id()] = false;
                            }
                        } else if (h.second > 0.9 * valid_read_length + read_infos[it->b_id()]->begin()) {
                            if (begin > h.first) {
                                is_valid_overlap[it->id()] = false;
                            }
                        }
                    }
                }
            }

            bool is_valid = it->update(
                read_infos[it->a_id()]->begin(),
                read_infos[it->a_id()]->end(),
                read_infos[it->b_id()]->begin(),
                read_infos[it->b_id()]->end()
            );

            if (!is_valid) {
                is_valid_overlap[it->id()] = false;
                continue;
            }

            it->set_type(classifyOverlap(it));
            switch (it->type()) {
                case 0:
                    is_valid_overlap[it->id()] = false;
                    break;
                case 1:
                    read_infos[it->a_id()].reset();
                    is_valid_read[it->a_id()] = false;
                    is_valid_overlap[it->id()] = false;
                    break;
                case 2:
                    read_infos[it->b_id()].reset();
                    is_valid_read[it->b_id()] = false;
                    is_valid_overlap[it->id()] = false;
                    break;
                default:
                    break;
            }
        }

        shrinkVector(overlaps, current_overlap_id, is_valid_overlap);

        if (status == false) {
            break;
        }
    }
    oreader.reset();

    // check if all non valid overlaps are deleted
    {
        bool shrink = false;
        for (const auto& it: overlaps) {
            if (!is_valid_read[it->a_id()] || !is_valid_read[it->b_id()]) {
                shrink = true;
                is_valid_overlap[it->id()] = false;
            }
        }
        if (shrink) {
            shrinkVector(overlaps, 0, is_valid_overlap);
        }
    }

    auto rreader = bioparser::createReader<Read, bioparser::FastqReader>(reads_path);
    while (true) {
        uint64_t current_read_id = reads.size();
        auto status = rreader->read_objects(reads, kChunkSize);

        for (uint64_t i = current_read_id; i < reads.size(); ++i) {
            auto& it = reads[i];
            if (it->id() >= is_valid_read.size()) {
                read_infos.resize(it->id() + 1);
                is_valid_read.resize(it->id() + 1, false);
                continue;
            }
            if (!is_valid_read[it->id()]) {
                continue;
            }
            it->trim_sequence(read_infos[it->id()]->begin(), read_infos[it->id()]->end());
        }

        shrinkVector(reads, current_read_id, is_valid_read);

        if (status == false) {
            break;
        }
    }

    uint64_t num_valid_reads = 0, num_valid_overlaps = 0;
    for (const auto& it: is_valid_read) if (it == true) ++num_valid_reads;
    for (const auto& it: is_valid_overlap) if (it == true) ++num_valid_overlaps;

    fprintf(stderr, "  number of valid overlaps = %lu (out of %lu)\n", num_valid_overlaps, is_valid_overlap.size());
    fprintf(stderr, "  number of valid reads = %lu (out of %lu)\n", num_valid_reads, is_valid_read.size());
    fprintf(stderr, "}\n\n");
}

void findChimericReads(const std::string& reads_path, const std::string& overlaps_path,
    uint32_t overlap_type) {

    fprintf(stderr, "Searching chimering reads (with reference) {\n");

    auto oreader = overlap_type == 0 ?
        bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
        bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

    std::vector<std::shared_ptr<Overlap>> overlaps;
    oreader->read_objects(overlaps, -1);

    auto rreader = bioparser::createReader<Read, bioparser::FastqReader>(reads_path);
    std::vector<std::shared_ptr<Read>> reads;
    rreader->read_objects(reads, -1);

    std::vector<bool> is_chimeric(reads.size(), false);

    uint32_t last_read_id = 0, last_overlap_id = 0, num_chimeras = 0, read_id = 1;
    float err = 0.25;

    for (const auto& it: overlaps) {
        if (last_read_id != it->a_id()) {
            uint32_t num_overlaps = it->id() - last_overlap_id;
            if (num_overlaps > 1) {

                // extend overlaps
                for (uint32_t i = last_overlap_id; i < it->id(); ++i) {
                    const auto& i_ovl = overlaps[i];

                    uint32_t i_begin = i_ovl->a_rc() ? i_ovl->a_length() - i_ovl->a_end() : i_ovl->a_begin();
                    uint32_t i_end = i_ovl->a_rc() ? i_ovl->a_length() - i_ovl->a_begin() : i_ovl->a_end();
                    uint32_t i_rbegin = i_ovl->b_rc() ? i_ovl->b_length() - i_ovl->b_end() : i_ovl->b_begin();
                    uint32_t i_rend = i_ovl->b_rc() ? i_ovl->b_length() - i_ovl->b_begin() : i_ovl->b_end();
                    uint32_t i_strand = i_ovl->a_rc() ^ i_ovl->b_rc();
                    uint32_t i_r = i_ovl->b_id();

                    while (true) {
                        bool is_changed = false;
                        for (uint32_t j = last_overlap_id; j < it->id(); ++j) {
                            if (j == i) {
                                continue;
                            }
                            const auto& j_ovl = overlaps[j];

                            uint32_t j_begin = j_ovl->a_rc() ? j_ovl->a_length() - j_ovl->a_end() : j_ovl->a_begin();
                            uint32_t j_end = j_ovl->a_rc() ? j_ovl->a_length() - j_ovl->a_begin() : j_ovl->a_end();

                            if (i_begin <= j_begin && j_end <= i_end) {
                                // contained repeat
                                continue;
                            }

                            uint32_t j_rbegin = j_ovl->b_rc() ? j_ovl->b_length() - j_ovl->b_end() : j_ovl->b_begin();
                            uint32_t j_rend = j_ovl->b_rc() ? j_ovl->b_length() - j_ovl->b_begin() : j_ovl->b_end();
                            uint32_t j_strand = j_ovl->a_rc() ^ j_ovl->b_rc();
                            uint32_t j_r = j_ovl->b_id();

                            int32_t diff_len = abs(j_begin > i_begin ? j_begin - i_end : i_begin - j_end) * (1 + err);
                            int32_t diff_rlen = abs(j_rbegin > i_rbegin ? j_rbegin - i_rend : i_rbegin - j_rend);

                            if (diff_rlen <= diff_len && i_strand == j_strand && i_r == j_r) {
                                // merge over low quality area
                                is_changed = true;

                                i_begin = std::min(i_begin, j_begin);
                                i_end = std::max(i_end, j_end);
                                i_rbegin = std::min(i_rbegin, j_rbegin);
                                i_rend = std::max(i_rend, j_rend);
                            }
                        }

                        if (!is_changed) {
                            break;
                        }

                        // update overlap
                        i_ovl->set_a_begin(i_ovl->a_rc() ? i_ovl->a_length() - i_end : i_begin);
                        i_ovl->set_a_end(i_ovl->a_rc() ? i_ovl->a_length() - i_begin : i_end);
                        i_ovl->set_b_begin(i_ovl->b_rc() ? i_ovl->b_length() - i_rend : i_rbegin);
                        i_ovl->set_b_end(i_ovl->b_rc() ? i_ovl->b_length() - i_rbegin : i_rend);
                    }
                }

                // check for repeats
                std::vector<bool> is_merged(num_overlaps, false);
                uint32_t num_merged_overlaps = 0;
                for (uint32_t i = last_overlap_id; i < it->id(); ++i) {
                    if (is_merged[i - last_overlap_id]) {
                        continue;
                    }

                    const auto& i_ovl = overlaps[i];
                    uint32_t i_begin = i_ovl->a_rc() ? i_ovl->a_length() - i_ovl->a_end() : i_ovl->a_begin();
                    uint32_t i_end = i_ovl->a_rc() ? i_ovl->a_length() - i_ovl->a_begin() : i_ovl->a_end();

                    for (uint32_t j = last_overlap_id; j < it->id(); ++j) {
                        if (j == i || is_merged[j - last_overlap_id]) {
                            continue;
                        }

                        const auto& j_ovl = overlaps[j];
                        uint32_t j_begin = j_ovl->a_rc() ? j_ovl->a_length() - j_ovl->a_end() : j_ovl->a_begin();
                        uint32_t j_end = j_ovl->a_rc() ? j_ovl->a_length() - j_ovl->a_begin() : j_ovl->a_end();

                        if (i_begin <= j_begin && j_end <= i_end) { // i_begin <= j_end && j_begin <= i_end) {
                            // contained/overlapped repeat
                            is_merged[j - last_overlap_id] = true;
                            ++num_merged_overlaps;
                            continue;
                        }
                    }
                }

                // check for circular genomes
                for (uint32_t i = last_overlap_id; i < it->id(); ++i) {
                    if (is_merged[i - last_overlap_id]) {
                        continue;
                    }

                    const auto& i_ovl = overlaps[i];

                    uint32_t i_begin = i_ovl->a_rc() ? i_ovl->a_length() - i_ovl->a_end() : i_ovl->a_begin();
                    uint32_t i_end = i_ovl->a_rc() ? i_ovl->a_length() - i_ovl->a_begin() : i_ovl->a_end();
                    uint32_t i_rbegin = i_ovl->b_rc() ? i_ovl->b_length() - i_ovl->b_end() : i_ovl->b_begin();
                    uint32_t i_rend = i_ovl->b_rc() ? i_ovl->b_length() - i_ovl->b_begin() : i_ovl->b_end();
                    uint32_t i_strand = i_ovl->a_rc() ^ i_ovl->b_rc();
                    uint32_t i_r = i_ovl->b_id();

                    for (uint32_t j = last_overlap_id; j < it->id(); ++j) {
                        if (j == i || is_merged[j - last_overlap_id]) {
                            continue;
                        }

                        const auto& j_ovl = overlaps[j];

                        uint32_t j_begin = j_ovl->a_rc() ? j_ovl->a_length() - j_ovl->a_end() : j_ovl->a_begin();
                        uint32_t j_end = j_ovl->a_rc() ? j_ovl->a_length() - j_ovl->a_begin() : j_ovl->a_end();
                        uint32_t j_rbegin = j_ovl->b_rc() ? j_ovl->b_length() - j_ovl->b_end() : j_ovl->b_begin();
                        uint32_t j_rend = j_ovl->b_rc() ? j_ovl->b_length() - j_ovl->b_begin() : j_ovl->b_end();
                        uint32_t j_strand = j_ovl->a_rc() ^ j_ovl->b_rc();
                        uint32_t j_r = j_ovl->b_id();

                        int32_t diff_len = abs(j_begin > i_begin ? j_begin - i_end : i_begin - j_end);
                        int32_t diff_rlen = abs(j_rbegin > i_rbegin ? it->b_length() - j_rend + i_rbegin : j_rbegin + it->b_length() - i_rend);
                        if (comparable(diff_rlen, diff_len, err) && i_strand == j_strand && i_r == j_r) {
                            // mergable
                            is_merged[j - last_overlap_id] = true;
                            ++num_merged_overlaps;
                        }
                    }
                }

                if (num_merged_overlaps + 1 != num_overlaps) {
                    // fprintf(stderr, "Chimeric: %u\n", last_read_id);
                    is_chimeric[last_read_id] = true;
                    ++num_chimeras;
                    for (uint32_t i = last_overlap_id; i < it->id(); ++i) {
                        if (is_merged[i - last_overlap_id]) {
                            continue;
                        }

                        const auto& i_ovl = overlaps[i];
                        uint32_t i_begin = i_ovl->a_rc() ? i_ovl->a_length() - i_ovl->a_end() : i_ovl->a_begin();
                        uint32_t i_end = i_ovl->a_rc() ? i_ovl->a_length() - i_ovl->a_begin() : i_ovl->a_end();

                        // fprintf(stderr, "Read id, begin, end = %d, %d, %d\n", last_read_id, i_begin, i_end);
                        fprintf(stdout, "@%u\n%s\n+\n%s\n", read_id++,
                            reads[last_read_id]->sequence().substr(i_begin, i_end - i_begin).c_str(),
                            reads[last_read_id]->quality().substr(i_begin, i_end - i_begin).c_str());
                    }
                }
            }

            last_read_id = it->a_id();
            last_overlap_id = it->id();
        }
    }

    for (uint32_t i = 0; i < reads.size(); ++i) {
        if (is_chimeric[i]) {
            continue;
        }
        fprintf(stdout, "@%u\n%s\n+\n%s\n", read_id++,
            reads[i]->sequence().c_str(),
            reads[i]->quality().c_str());
    }

    fprintf(stderr, "  found %u chimeric reads\n", num_chimeras);
    fprintf(stderr, "}\n");
}

void findUncontainedReads(const std::string& reads_path, const std::string& overlaps_path,
    uint32_t overlap_type) {

    fprintf(stderr, "Searching uncontained reads {\n");

    std::vector<std::shared_ptr<Read>> reads;
    auto rreader = bioparser::createReader<Read, bioparser::FastqReader>(reads_path);

    std::vector<uint32_t> read_lengths;

    while (true) {
        auto status = rreader->read_objects(reads, kChunkSize);
        read_lengths.resize(read_lengths.size() + reads.size());

        for (const auto& it: reads) {
            read_lengths[it->id()] = it->sequence().size() * kMinMappingPerc;
        }
        reads.clear();

        if (status == false) {
            break;
        }
    }
    rreader.reset();

    std::vector<uint32_t> read_type(read_lengths.size(), 0);

    std::vector<std::shared_ptr<Overlap>> overlaps;
    auto oreader = overlap_type == 0 ?
        bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
        bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

    uint32_t internal = 0, cont = 0, prefsuf = 0;
    while (true) {
        auto status = oreader->read_objects(overlaps, kChunkSize);

        for (const auto& it: overlaps) {
            it->set_type(classifyOverlap(it));
            read_type[it->a_id()] = std::max(read_type[it->a_id()], it->type() + 1);
            switch (it->type()) {
                case 0:
                    ++internal;
                    break;
                case 1:
                case 2:
                    ++cont;
                    break;
                case 3:
                case 4:
                default:
                    ++prefsuf;
                    break;
            }
        }

        overlaps.clear();

        if (status == false) {
            break;
        }
    }

    fprintf(stderr, "Int = %d, Con = %d, PS = %d\n", internal, cont, prefsuf);
    uint32_t num_uncontained_reads = 0;

    rreader = bioparser::createReader<Read, bioparser::FastqReader>(reads_path);
    while (true) {
        auto status = rreader->read_objects(reads, kChunkSize);

        for (const auto& it: reads) {
            if (read_type[it->id()] == 0 || read_type[it->id()] > 3) {
                fprintf(stdout, "@%s\n%s\n+\n%s\n",
                    it->name().c_str(),
                    it->sequence().c_str(),
                    it->quality().c_str());
                ++num_uncontained_reads;
            }
        }
        reads.clear();

        if (status == false) {
            break;
        }
    }

    fprintf(stderr, "  found %u uncontained reads (out of %lu)\n", num_uncontained_reads, read_type.size());
    fprintf(stderr, "}\n");
}

void fastaToFastq(const std::string& reads_path) {

    std::vector<std::shared_ptr<Read>> reads;
    auto rreader = bioparser::createReader<Read, bioparser::FastaReader>(reads_path);

    std::vector<uint32_t> read_lengths;

    while (true) {
        auto status = rreader->read_objects(reads, kChunkSize);
        read_lengths.resize(read_lengths.size() + reads.size());

        for (const auto& it: reads) {
            fprintf(stdout, "@%s\n%s\n+\n%s\n",
                it->name().c_str(),
                it->sequence().c_str(),
                std::string(it->sequence().size(), '*').c_str());
        }
        reads.clear();

        if (status == false) {
            break;
        }
    }
}

bool comparable(double a, double b, double eps) {
    return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
        (b >= a * (1 - eps) && b <= a * (1 + eps));
}

}
