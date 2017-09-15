/*!
 * @file preprocess.cpp
 *
 * @brief Preprocessing source file
 */
/*
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <set>

#include "read.hpp"
#include "overlap.hpp"
#include "timer.hpp"
#include "preprocess.hpp"

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"

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
            if (is_valid_overlap[overlaps[j]->id()] == false ||
                overlaps[i]->b_id() != overlaps[j]->b_id()) {
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
    if (overlap->a_begin() <= overlap->b_begin() &&
        (overlap->a_length() - overlap->a_end()) <= (overlap->b_length() - overlap->b_end())) {
        return 1; // a contained
    }
    if (overlap->a_begin() >= overlap->b_begin() &&
        (overlap->a_length() - overlap->a_end()) >= (overlap->b_length() - overlap->b_end())) {
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

    auto save_valid_overlaps = [&](void) -> void {
        num_duplicate_overlaps += removeDuplicateOverlaps(is_valid_overlap, current_overlaps);

        // save valid overlaps
        for (uint32_t i = 0; i < current_overlaps.size(); ++i) {
            if (is_valid_overlap[current_overlaps[i]->id()]) {
                const auto& o = current_overlaps[i];

                mappings[o->a_id()].push_back(((o->a_rc() ? o->a_length() - o->a_end() : o->a_begin()) + 1) << 1 | 0);
                mappings[o->a_id()].push_back(((o->a_rc() ? o->a_length() - o->a_begin() : o->a_end()) - 1) << 1 | 1);

                mappings[o->b_id()].push_back(((o->b_rc() ? o->b_length() - o->b_end() : o->b_begin()) + 1) << 1 | 0);
                mappings[o->b_id()].push_back(((o->b_rc() ? o->b_length() - o->b_begin() : o->b_end()) - 1) << 1 | 1);
            }
        }

        current_overlaps.clear();
    };

    fprintf(stderr, "  Loading overlaps {\n");
    Timer load_timer;
    load_timer.start();

    while (true) {
        auto status = oreader->read_objects(overlaps, kChunkSize);
        is_valid_overlap.resize(is_valid_overlap.size() + overlaps.size(), true);

        for (const auto& it: overlaps) {

            uint32_t max_read_id = std::max(it->a_id(), it->b_id());
            if (mappings.size() <= max_read_id) {
                mappings.resize(max_read_id + 1);
                read_lengths.resize(max_read_id + 1);
            }

            if (read_infos.size() <= max_read_id) {
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

            if (current_overlaps.size() != 0 && current_overlaps.front()->a_id() != it->a_id()) {
                save_valid_overlaps();
            }
            current_overlaps.push_back(it);
        }

        overlaps.clear();

        if (status == false) {
            save_valid_overlaps();
        }

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

    load_timer.stop();
    fprintf(stderr, "    number of self overlaps = %u\n", num_self_overlaps);
    fprintf(stderr, "    number of duplicate overlaps = %u\n", num_duplicate_overlaps);
    load_timer.print("    time =");
    fprintf(stderr, "  }\n");

    uint32_t num_prefiltered_reads = 0;
    std::vector<bool> is_valid_read(read_infos.size(), true);
    for (uint32_t i = 0; i < read_infos.size(); ++i) {
        if (read_infos[i] == nullptr) {
            is_valid_read[i] = false;
            ++num_prefiltered_reads;
        }
    }

    // find coverage medians
    {
        fprintf(stderr, "  Calculating coverage madians {\n");
        Timer timer;
        timer.start();

        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) {
                    read_infos[id]->find_coverage_median();
                }, i));
        }
        for (const auto& it: thread_futures) {
            it.wait();
        }

        std::vector<uint16_t> medians;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }
            medians.emplace_back(read_infos[i]->coverage_median());
        }

        std::sort(medians.begin(), medians.end());
        dataset_coverage_median = medians.size() % 2 == 1 ? medians[medians.size() / 2] :
            (medians[medians.size() / 2 - 1] + medians[medians.size() / 2]) / 2;

        timer.stop();
        fprintf(stderr, "    median of coverage medians = %u\n",
            (uint32_t) dataset_coverage_median);
        timer.print("    time =");
        fprintf(stderr, "  }\n");
    }

    // find chimeric reads by looking for tight coverage pits
    {
        fprintf(stderr, "  Chimeric check {\n");
        Timer timer;
        timer.start();

        std::vector<std::future<bool>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) {
                    return read_infos[id]->find_coverage_pits(1.817,
                        kSlopeWidth, kSlopeWidthRatio, dataset_coverage_median);
                }, i));
        }
        uint32_t num_chimeric_reads = 0;
        for (auto& it: thread_futures) {
            it.wait();
            if (it.get() == true) {
                ++num_chimeric_reads;
            }
        }

        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }
            if (!read_infos[i]->is_valid()) {
                read_infos[i].reset();
                is_valid_read[i] = false;
            }
        }

        timer.stop();
        fprintf(stderr, "    number of chimeric reads = %u\n", num_chimeric_reads);
        timer.print("    time:");
        fprintf(stderr, "  }\n");
    }

    // find longest contiguous read region which has coverage larger than kMinCoverage
    {
        fprintf(stderr, "  Trimming data {\n");
        Timer timer;
        timer.start();

        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                is_valid_read[i] = false;
                ++num_prefiltered_reads;
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) {
                    read_infos[id]->find_valid_region(kMinCoverage);
                }, i));
        }
        for (auto& it: thread_futures) {
            it.wait();
        }

        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            if (!read_infos[i]->is_valid()) {
                read_infos[i].reset();
                is_valid_read[i] = false;
                ++num_prefiltered_reads;
            }
        }

        timer.stop();
        fprintf(stderr, "    number of filtered reads = %u\n", num_prefiltered_reads);
        timer.print("    time =");
        fprintf(stderr, "  }\n");
    }

    // filter low median reads
    {
        fprintf(stderr, "  Filtering low median reads {\n");
        Timer timer;
        timer.start();

        uint32_t num_low_median_reads = 0;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            if (read_infos[i]->coverage_median() * kMedianRatio < dataset_coverage_median) {
                read_infos[i].reset();
                is_valid_read[i] = false;
                ++num_low_median_reads;
            }
        }

        timer.stop();
        fprintf(stderr, "    number of reads with coverage median < %u = %u\n",
            (uint32_t) (dataset_coverage_median / kMedianRatio), num_low_median_reads);
        timer.print("    time =");
        fprintf(stderr, "  }\n");
    }

    // correct coverage graphs
    {
        fprintf(stderr, "  Correcting coverage graphs {\n");
        Timer timer;
        timer.start();

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

        // update coverage medians
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
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr) {
                continue;
            }

            medians.emplace_back(read_infos[i]->coverage_median());
        }

        std::sort(medians.begin(), medians.end());
        dataset_coverage_median = medians.size() % 2 == 1 ? medians[medians.size() / 2] :
            (medians[medians.size() / 2 - 1] + medians[medians.size() / 2]) / 2;

        timer.stop();
        fprintf(stderr, "    updated median of coverage medians = %u\n",
            (uint32_t) dataset_coverage_median);
        timer.print("    time =");
        fprintf(stderr, "  }\n");
    }

    // smooth coverage graphs
    {
        fprintf(stderr, "  Smoothing coverage graphs {\n");
        Timer timer;
        timer.start();

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

        timer.stop();
        timer.print("    time =");
        fprintf(stderr, "  }\n");
    }

    // find coverage hills which represent repetitive regions
    {
        fprintf(stderr, "  Hill detection {\n");
        Timer timer;
        timer.start();

        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < read_infos.size(); ++i) {
            if (read_infos[i] == nullptr || !read_infos[i]->coverage_hills().empty()) {
                continue;
            }

            thread_futures.emplace_back(thread_pool->submit_task(
                [&](uint32_t id) {
                    read_infos[id]->find_coverage_hills(kSlopeRatio, kSlopeWidth,
                        kSlopeWidthRatio, kHillWidthRatio, dataset_coverage_median);
                }, i));
        }
        for (auto& it: thread_futures) {
            it.wait();
        }

        uint32_t num_left_repeats = 0, num_right_repeats = 0;
        for (const auto& it: read_infos) {
            if (it == nullptr) {
                continue;
            }

            uint32_t valid_read_length = it->end() - it->begin();
            for (const auto& hill: it->coverage_hills()) {
                if (hill.first < 0.1 * valid_read_length + it->begin()) {
                    ++num_left_repeats;
                } else if (hill.second > 0.9 * valid_read_length + it->begin()) {
                    ++num_right_repeats;
                }
            }
        }

        timer.stop();
        fprintf(stderr, "    number of left repeat reads = %u\n", num_left_repeats);
        fprintf(stderr, "    number of right repeat reads = %u\n", num_right_repeats);
        timer.print("    time =");
        fprintf(stderr, "  }\n");
    }

    std::vector<uint32_t> specials = {};
    for (const auto& id: specials) {
        if (read_infos[id] == nullptr) {
            continue;
        }
        read_infos[id]->print_csv("graphs/h" + std::to_string(id), dataset_coverage_median);
    }

    load_timer.reset();
    load_timer.start();
    fprintf(stderr, "  Loading reads and overlaps {\n");

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
                uint32_t fuzz = 0.05 * valid_read_length;

                for (const auto& h: read_infos[it->a_id()]->coverage_hills()) {
                    if (begin < h.second && h.first < end) {
                        if (h.first < 0.1 * valid_read_length + read_infos[it->a_id()]->begin()) {
                            if (end < h.second + 0.05 * valid_read_length) {
                                is_valid_overlap[it->id()] = false;
                            }
                        } else if (h.second > 0.9 * valid_read_length + read_infos[it->a_id()]->begin()) {
                            if (begin > h.first - 0.05 * valid_read_length) {
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
                uint32_t fuzz = 0.05 * valid_read_length;

                for (const auto& h: read_infos[it->b_id()]->coverage_hills()) {
                    if (begin < h.second && h.first < end) {
                        if (h.first < 0.1 * valid_read_length + read_infos[it->b_id()]->begin()) {
                            if (end < h.second + fuzz) {
                                is_valid_overlap[it->id()] = false;
                            }
                        } else if (h.second > 0.9 * valid_read_length + read_infos[it->b_id()]->begin()) {
                            if (begin > h.first - fuzz) {
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

    auto rreader = bioparser::createReader<Read, bioparser::FastaReader>(reads_path);
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

    load_timer.stop();
    fprintf(stderr, "    number of valid overlaps = %lu (out of %lu)\n", num_valid_overlaps, is_valid_overlap.size());
    fprintf(stderr, "    number of valid reads = %lu (out of %lu)\n", num_valid_reads, is_valid_read.size());
    load_timer.print("    time =");
    fprintf(stderr, "  }\n}\n\n");
}

}
*/
