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

// prefilter
constexpr uint32_t kMinCoverage = 3;

// preprocess
constexpr double kMaxOverhangToOverlapRatio = 0.875;
constexpr double kSlopeRatio = 1.3;
constexpr uint32_t kSlopeWidth = 500;
constexpr double kSlopeWidthRatio = 0.05;
constexpr double kHillWidthRatio = 0.7;
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
                is_valid_overlap[overlaps[i]->id()] = false;
                break;
            } else {
                is_valid_overlap[overlaps[j]->id()] = false;
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

// call before mappings sort!
void printPileGraph(const std::vector<uint32_t>& mappings, const std::string& path) {

    std::ofstream out(path);
    for (uint32_t i = 0; i < mappings.size(); i+= 2) {
        out << (mappings[i] >> 1) << " " << (mappings[i + 1] >> 1) << std::endl;
    }
    out.close();
};

void printCoverageGraph(std::vector<uint32_t>& mappings, const std::string& path, double dataset_median) {

    if (mappings.empty()) {
        return;
    }
    std::sort(mappings.begin(), mappings.end());

    std::vector<uint32_t> coverage_graph((mappings.back() >> 1) + 1, 0);
    int32_t coverage = 0;
    uint32_t last_m = 0, total = 0;
    for (const auto& m: mappings) {
        if (coverage > 0) {
            for (uint32_t i = last_m; i < (m >> 1); ++i) {
                coverage_graph[i] += coverage;
                total += coverage;
            }
        }
        last_m = m >> 1;
        if (m & 1){
            --coverage;
        } else {
            ++coverage;
        }
    }

    std::vector<uint32_t> tmp(coverage_graph.begin() + (mappings.front()>>1), coverage_graph.end());
    std::sort(tmp.begin(), tmp.end());
    double median = tmp.size() % 2 == 1 ? tmp[tmp.size() / 2] : (tmp[tmp.size() / 2 - 1] + tmp[tmp.size() / 2]) / (double) 2;

    std::ofstream out(path);
    for (uint32_t i = 0; i < coverage_graph.size(); ++i) {
        out << i << " " << coverage_graph[i] << " " << 0 << " " << median << " " << dataset_median << std::endl;
    }
    out.close();
}

bool findRegions(std::vector<std::pair<uint32_t, uint32_t>>& valid_regions,
    std::pair<uint32_t, uint32_t>& longest_region, std::vector<uint32_t>& mappings) {

    if (mappings.empty()) {
        return false;
    }
    std::sort(mappings.begin(), mappings.end());

    int32_t coverage = 0, min_coverage = kMinCoverage, region_begin = 0;
    for (const auto& m: mappings) {
        int32_t old_coverage = coverage;
        if (m & 1) {
            --coverage;
        } else {
            ++coverage;
        }
        int32_t pos = m >> 1;
        if (old_coverage < min_coverage && coverage >= min_coverage) {
            region_begin = pos;
        } else if (old_coverage >= min_coverage && coverage < min_coverage) {
            if (pos - region_begin > 500) {
                valid_regions.emplace_back(region_begin, pos);
            }
        }
    }

    if (valid_regions.empty()) {
        return false;
    } else {
        longest_region = valid_regions.front();
        for (uint32_t i = 1; i < valid_regions.size(); ++i) {
            if (valid_regions[i].second - valid_regions[i].first > longest_region.second - longest_region.first) {
                longest_region = valid_regions[i];
            }
        }
    }
    return true;
}

void findHills(std::vector<std::pair<uint32_t, uint32_t>>& hills,
    std::vector<uint32_t>& mappings, double dataset_median, uint32_t id, bool print) {

    if (mappings.empty()) {
        return;
    }
    std::sort(mappings.begin(), mappings.end());

    std::vector<uint32_t> coverage_graph((mappings.back() >> 1) + 1, 0);
    int32_t coverage = 0;
    uint32_t last_m = 0;
    for (const auto& m: mappings) {
        if (coverage > 0) {
            for (uint32_t i = last_m; i < (m >> 1); ++i) {
                coverage_graph[i] += coverage;
            }
        }
        last_m = m >> 1;
        if (m & 1){
            --coverage;
        } else {
            ++coverage;
        }
    }

    std::vector<uint32_t> tmp(coverage_graph.begin() + (mappings.front()>>1), coverage_graph.end());
    std::sort(tmp.begin(), tmp.end());
    double median = tmp.size() % 2 == 1 ? tmp[tmp.size() / 2] : (tmp[tmp.size() / 2 - 1] + tmp[tmp.size() / 2]) / (double) 2;

    auto window_add = [](std::deque<std::pair<int32_t, int32_t>>& window,
        int32_t value, int32_t position) -> void {

        while (!window.empty() && window.back().second <= value) {
            window.pop_back();
        }
        window.emplace_back(position, value);
        return;
    };

    auto window_update = [](std::deque<std::pair<int32_t, int32_t>>& window,
        int32_t position) -> void {

        while (!window.empty() && window.front().first <= position) {
            window.pop_front();
        }
        return;
    };

    std::deque<std::pair<int32_t, int32_t>> left_window, right_window;
    std::vector<int32_t> slopes;

    int32_t k = std::max(kSlopeWidth, uint32_t(kSlopeWidthRatio * ((mappings.back() >> 1) - (mappings.front() >> 1))));
    int32_t read_length = coverage_graph.size();

    uint32_t beg = mappings.front()>>1, len = (mappings.back()>>1) - (mappings.front()>>1);

    for (int32_t i = -1 * k + 2; i < read_length; ++i) {
        if (i < read_length - k) {
            window_add(right_window, coverage_graph[i + k], i + k);
        }
        window_update(right_window, i);

        if (i > 0) {
            window_add(left_window, coverage_graph[i - 1], i - 1);
            window_update(left_window, i - 1 - k);

            // if (coverage_graph[i - 1] < dataset_median && (i == read_length - 1 || coverage_graph[i + 1] < dataset_median)) continue;
            int32_t current = coverage_graph[i] * kSlopeRatio;
            if (coverage_graph[i - 1] > dataset_median && left_window.front().second > dataset_median && left_window.front().second > current) {
                slopes.push_back(i << 1 | 0);
            }
            if ((i == read_length - 1 || coverage_graph[i + 1] > dataset_median) && !right_window.empty() && right_window.front().second > dataset_median && right_window.front().second > current) {
                slopes.push_back(i << 1 | 1);
            }
        }
    }

    uint32_t max_width = ((mappings.back() >> 1) - (mappings.front() >> 1)) * kHillWidthRatio;

    if (slopes.size() > 2) {

        uint32_t ldownslope = 0, fdownslope = 0;
        bool found_fds = false;

        uint32_t fupslope = 0, lupslope = 0;
        bool found_fus = false;

        std::vector<std::pair<uint32_t, uint32_t>> slope_regions;

        uint32_t slope_width = k; // 2 * k;
        for (uint32_t s = 0; s < slopes.size(); ++s) {
            if (slopes[s] & 1) {
                if (found_fus) {
                    if ((slopes[s] >> 1) - fupslope > slope_width) {
                        slope_regions.emplace_back(fupslope << 1 | 1, lupslope);
                        fupslope = slopes[s] >> 1;
                        lupslope = fupslope;
                    } else {
                        lupslope = slopes[s] >> 1;
                    }
                } else {
                    found_fus = true;
                    fupslope = slopes[s] >> 1;
                    lupslope = fupslope;
                }
            } else {
                if (found_fds) {
                    if ((slopes[s] >> 1) - fdownslope > slope_width) {
                        slope_regions.emplace_back(fdownslope << 1 | 0, ldownslope);
                        fdownslope = slopes[s] >> 1;
                        ldownslope = fdownslope;
                    } else {
                        ldownslope = slopes[s] >> 1;
                    }
                } else {
                    found_fds = true;
                    fdownslope = slopes[s] >> 1;
                    ldownslope = fdownslope;
                }
            }
        }

        if (found_fus) slope_regions.emplace_back(fupslope << 1 | 1, lupslope);
        if (found_fds) slope_regions.emplace_back(fdownslope << 1 | 0, ldownslope);
        std::sort(slope_regions.begin(), slope_regions.end());

        // for (const auto& it: slope_regions) {
        //    fprintf(stderr, "%d (%d %d) \n", it.first & 1, it.first >> 1, it.second);
        // }

        // chimeric check - TODO: rerranging slopes does not work all the time
        auto rearrange_slopes = [&](uint32_t i, uint32_t j) -> void {
            uint32_t begin = std::max(slope_regions[i].first >> 1, slope_regions[j].first >> 1);
            uint32_t end = std::min(slope_regions[i].second, slope_regions[j].second);
            uint32_t min_left_id = begin, min_right_id = begin;
            for (uint32_t s = begin + 1; s < end; ++s) {
                if (coverage_graph[s] < coverage_graph[min_left_id]) {
                    min_left_id = s;
                }
                if (coverage_graph[s] <= coverage_graph[min_right_id]) {
                    min_right_id = s;
                }
            }
            ++min_left_id;
            --min_right_id;

            slope_regions[i].first = min_left_id << 1 | 0;
            slope_regions[i].second = min_left_id;
            slope_regions[j].first = min_right_id << 1 | 1;
            slope_regions[j].second = min_right_id;
        };

        for (uint32_t s = 0; s < slope_regions.size() - 1; ++s) {
            if ((slope_regions[s].first & 1) && !(slope_regions[s + 1].first & 1) &&
                slope_regions[s].second > (slope_regions[s + 1].first >> 1)) {
                rearrange_slopes(s, s + 1);
            }
        }

        for (uint32_t s = 0; s < slope_regions.size() - 1; ++s) {
            if (slope_regions[s].second > (slope_regions[s + 1].first >> 1)) {
                rearrange_slopes(s, s + 1);
            }
        }
        // fprintf(stderr, "Chimeric check\n");
        /* for (const auto& it: slope_regions) {
            fprintf(stderr, "%d (%d %d) \n", it.first & 1, it.first >> 1, it.second);
        }
        fprintf(stderr, "\n"); */

        std::vector<uint32_t> slope_graph(coverage_graph.size(), 0);
        slopes.clear();

        // find hills
        auto check_hill = [&](uint32_t begin, uint32_t end, double median) -> bool {
            for (uint32_t i = begin; i < end; ++i) {
                if (coverage_graph[i] < median) {
                    return false;
                }
            }
            return true;
        };

        for (uint32_t r = 0; r < slope_regions.size() - 1; ++r) {
            if ((slope_regions[r].first & 1) && !(slope_regions[r + 1].first & 1)) {
                if (slope_regions[r + 1].second - (slope_regions[r].first >> 1) < max_width) {

                    slope_graph[(slope_regions[r].first >> 1)] = 2;
                    slope_graph[(slope_regions[r + 1].second)] = 1;
                    // print = true;

                    // fprintf(stderr, "%d: (%d,%d) - (%d,%d)\n", id, slope_regions[r].first >> 1, slope_regions[r].second,
                    //    slope_regions[r + 1].first >> 1, slope_regions[r + 1].second);

                    if (check_hill(slope_regions[r].second + 1, (slope_regions[r + 1].first >> 1) - 1, dataset_median) == false) {
                        continue;
                    }

                    if ((slope_regions[r].first >> 1) < 0.075 * len + beg) {
                        // left hill
                        hills.emplace_back(slope_regions[r].first >> 1, slope_regions[r + 1].second);
                    } else if (slope_regions[r + 1].second > 0.925 * len + beg) {
                        // right hill
                        hills.emplace_back(slope_regions[r].first >> 1, slope_regions[r + 1].second);
                    }
                }
            }
        }

        std::set<uint32_t> special = {}; // {30315, 25972, 12512, 43069, 28793, 27602, 8291, 40019, 22178, 18285, 30224, 30669, 42382, 28437, 2579};

        if (print || special.count(id) != 0) {
            std::ofstream out("graphs/h" + std::to_string(id));
            for (uint32_t i = 0; i < coverage_graph.size(); ++i) {
                out << i << " " << coverage_graph[i] << " " << slope_graph[i] << " " << median << " " << dataset_median << std::endl;
            }
            out.close();
        }
    }
}

void findHillsWithReference(std::vector<std::pair<uint32_t, uint32_t>>& hills,
    std::vector<uint32_t>& mappings, uint32_t id, bool print) {

    if (mappings.empty()) {
        return;
    }
    std::sort(mappings.begin(), mappings.end());

    std::vector<uint32_t> coverage_graph((mappings.back() >> 1) + 1, 0);
    int32_t coverage = 0, min_coverage = 2, region_begin = 0;
    uint32_t last_m = 0;
    for (const auto& m: mappings) {
        if (coverage > 0) {
            for (uint32_t i = last_m; i < (m >> 1); ++i) {
                coverage_graph[i] += coverage;
            }
        }
        last_m = m >> 1;
        int32_t old_coverage = coverage;
        if (m & 1){
            --coverage;
        } else {
            ++coverage;
        }
        int32_t pos = m >> 1;
        if (old_coverage < min_coverage && coverage >= min_coverage) {
            region_begin = pos;
        } else if (old_coverage >= min_coverage && coverage < min_coverage) {
            hills.emplace_back(region_begin, pos);
        }
    }

    if (!hills.empty() && print) {
        std::vector<uint32_t> slope_graph(coverage_graph.size(), 0);
        for (const auto& it: hills) {
            slope_graph[it.first] = 2;
            slope_graph[it.second] = 1;
        }

        std::ofstream out("graphs/c" + std::to_string(id));
        for (uint32_t i = 0; i < coverage_graph.size(); ++i) {
            out << i << " " << coverage_graph[i] << " " << slope_graph[i] << std::endl;
        }
        out.close();
    }
}

bool findPits(std::vector<uint32_t>& mappings, uint32_t id, bool print) {

    if (mappings.empty()) {
        return false;
    }
    std::sort(mappings.begin(), mappings.end());

    auto deque_add = [](std::deque<std::pair<int32_t, int32_t>>& window,
        int32_t value, int32_t position) -> void {

        while (!window.empty() && window.back().second <= value) {
            window.pop_back();
        }
        window.emplace_back(position, value);
        return;
    };

    auto deque_update = [](std::deque<std::pair<int32_t, int32_t>>& window,
        int32_t position) -> void {

        while (!window.empty() && window.front().first <= position) {
            window.pop_front();
        }
        return;
    };

    std::vector<uint32_t> coverage_graph((mappings.back() >> 1) + 1, 0);
    int32_t coverage = 0;
    uint32_t last_m = 0;
    for (const auto& m: mappings) {
        if (coverage > 0) {
            for (uint32_t i = last_m; i < (m >> 1); ++i) {
                coverage_graph[i] += coverage;
            }
        }
        last_m = m >> 1;
        if (m & 1) {
            --coverage;
        } else {
            ++coverage;
        }
    }

    std::deque<std::pair<int32_t, int32_t>> left_window, right_window;
    std::vector<int32_t> slopes;

    int32_t k = std::max(kSlopeWidth, uint32_t(kSlopeWidthRatio *
        ((mappings.back() >> 1) - (mappings.front() >> 1))));
    int32_t length = coverage_graph.size();

    for (int32_t i = -1 * k + 2; i < length - 1; ++i) {
        if (i < length - k) {
            deque_add(right_window, coverage_graph[i + k], i + k);
        }
        deque_update(right_window, i);

        if (i > 0) {
            deque_add(left_window, coverage_graph[i - 1], i - 1);
            deque_update(left_window, i - 1 - k);

            int32_t current = coverage_graph[i] * 1.817;
            if (left_window.front().second > current) {
                slopes.push_back(i << 1 | 0);
            }
            if (!right_window.empty() && right_window.front().second > current) {
                slopes.push_back(i << 1 | 1);
            }
        }
    }

    bool is_chimeric = false;
    if (slopes.size() > 0) {
        std::sort(slopes.begin(), slopes.end());
        for (uint32_t i = 0; i < slopes.size() - 1; ++i) {
            if (!(slopes[i] & 1) && (slopes[i + 1] & 1) && (slopes[i + 1] >> 1) - (slopes[i] >> 1) < k) {
                is_chimeric = true;
                break;
            }
        }
    }

    if (print && is_chimeric) {
        std::ofstream out("graphs/c" + std::to_string(id));
        for (uint32_t i = 0; i < coverage_graph.size(); ++i) {
            out << i << " " << coverage_graph[i] << " " <<  0 << std::endl;
        }
        out.close();
    }

    return is_chimeric;
}

void preprocessData(std::vector<std::shared_ptr<Read>>& reads, std::vector<std::shared_ptr<Overlap>>& overlaps,
    const std::string& reads_path, const std::string& overlaps_path, const std::string& mappings_path,
    uint32_t overlap_type, std::shared_ptr<thread_pool::ThreadPool> thread_pool, bool prefilter) {

    fprintf(stderr, "Preprocessing data {\n");

    std::vector<std::shared_ptr<Overlap>> current_overlaps;
    std::vector<std::vector<uint32_t>> mappings;
    std::vector<bool> is_valid_overlap;
    uint32_t num_duplicate_overlaps = 0;
    uint32_t num_self_overlaps = 0;

    std::vector<bool> is_valid_read;

    auto oreader = overlap_type == 0 ?
        bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
        bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

    while (true) {
        auto status = oreader->read_objects(overlaps, kChunkSize);
        is_valid_overlap.resize(is_valid_overlap.size() + overlaps.size(), true);

        for (const auto& it: overlaps) {

            uint32_t max_read_id = std::max(it->a_id(), it->b_id());
            if (mappings.size() <= max_read_id) {
                mappings.resize(max_read_id + 1);
            }

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

        if (status == false) {
            break;
        }
    }
    oreader.reset();
    is_valid_read.resize(mappings.size(), true);

    fprintf(stderr, "  number of self overlaps = %u\n", num_self_overlaps);
    fprintf(stderr, "  number of duplicate overlaps = %u\n", num_duplicate_overlaps);

    // find contiguous read regions which have coverage larger than kMinCoverage
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> valid_regions;
    std::vector<std::pair<uint32_t, uint32_t>> longest_region;
    if (prefilter) {
        fprintf(stderr, "  Prefiltering data {\n");

        valid_regions.resize(mappings.size());
        longest_region.resize(mappings.size());

        std::vector<std::future<bool>> thread_futures;
        for (uint32_t i = 0; i < mappings.size(); ++i) {
            thread_futures.emplace_back(thread_pool->submit_task(findRegions,
                std::ref(valid_regions[i]), std::ref(longest_region[i]),
                std::ref(mappings[i])));
        }
        for (uint32_t i = 0; i < mappings.size(); ++i) {
            thread_futures[i].wait();
            is_valid_read[i] = thread_futures[i].get();
        }

        uint32_t num_prefiltered_reads = 0;
        for (const auto& it: is_valid_read) if (it == false) ++num_prefiltered_reads;

        fprintf(stderr, "    number of prefiltered reads = %u\n", num_prefiltered_reads);
        fprintf(stderr, "  }\n");
    }

    std::vector<double> medians;
    uint64_t median_sum = 0;
    std::vector<double> read_medians(mappings.size(), 0);
    // for (const auto& it: mappings) {
    for (uint32_t i = 0; i < mappings.size(); ++i) {
        const auto& it = mappings[i];
        if (it.empty()) continue;
        std::vector<uint32_t> coverage_graph((it.back() >> 1) + 1, 0);
        int32_t coverage = 0;
        uint32_t last_m = 0;
        for (const auto& m: it) {
            if (coverage > 0) {
                for (uint32_t i = last_m; i < (m >> 1); ++i) {
                    coverage_graph[i] += coverage;
                }
            }
            last_m = m >> 1;
            if (m & 1){
                --coverage;
            } else {
                ++coverage;
            }
        }

        std::vector<uint32_t> cov(coverage_graph.begin()+(it.front()>>1), coverage_graph.end());
        if (cov.empty()) {
            continue;
        }

        std::sort(cov.begin(), cov.end());
        double median = cov.size() % 2 == 1 ? cov[cov.size() / 2] : (cov[cov.size() / 2 - 1] + cov[cov.size() / 2]) / (double) 2;
        medians.push_back(median);
        median_sum += median;
        read_medians[i] = median;
    }

    std::sort(medians.begin(), medians.end());
    double medians_median = (medians.size() % 2 == 1 ? medians[medians.size() / 2] : (medians[medians.size() / 2 - 1] + medians[medians.size() / 2]) / (double) 2);
    fprintf(stderr, "Median of medians %lf\n", medians_median);
    double average_median = median_sum / (double) medians.size();
    fprintf(stderr, "Avg median = %lf\n", average_median);
    medians.clear();

    uint32_t bad_reads = 0;
    for (uint32_t i = 0; i < mappings.size(); ++i) {
        if (read_medians[i] * kMedianRatio < medians_median) {
            std::vector<uint32_t>().swap(mappings[i]);
            is_valid_read[i] = false;
            ++bad_reads;
        }
    }
    fprintf(stderr, "Reads with median < median of medians / %lf = %u (out of %zu)\n", kMedianRatio, bad_reads, mappings.size());

    std::vector<uint32_t> specials = {};
    for (const auto& it: specials) {
        printCoverageGraph(mappings[it], "graphs/h" + std::to_string(it), average_median);
    }

    // find chimeric reads
    {
        std::vector<std::future<bool>> thread_futures;
        for (uint32_t i = 0; i < mappings.size(); ++i) {
            thread_futures.emplace_back(thread_pool->submit_task(findPits,
                std::ref(mappings[i]), i, false));
        }
        for (uint32_t i = 0; i < mappings.size(); ++i) {
            thread_futures[i].wait();
            if (thread_futures[i].get() == true) {
                std::vector<uint32_t>().swap(mappings[i]);
                is_valid_read[i] = false;
            }
        }
    }

    // find coverage hills which represent repeats or are parts of chimeric reads
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> hill_regions;
    std::vector<uint8_t> read_types(mappings.size());
    if (mappings_path.empty()) {
        hill_regions.resize(mappings.size());
        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < mappings.size(); ++i) {
            thread_futures.emplace_back(thread_pool->submit_task(findHills,
                std::ref(hill_regions[i]), std::ref(mappings[i]),
                average_median, i, false));
        }
        for (uint32_t i = 0; i < mappings.size(); ++i) {
            thread_futures[i].wait();
            std::vector<uint32_t>().swap(mappings[i]);
        }
        std::vector<std::vector<uint32_t>>().swap(mappings);
    } else {
        oreader = overlap_type == 0 ?
            bioparser::createReader<Overlap, bioparser::MhapReader>(mappings_path) :
            bioparser::createReader<Overlap, bioparser::PafReader>(mappings_path);
        oreader->read_objects(overlaps, -1);

        std::vector<std::vector<uint32_t>> ref_mappings;
        for (const auto& it: overlaps) {
            if (ref_mappings.size() <= it->a_id()) {
                ref_mappings.resize(it->a_id() + 1);
            }

            ref_mappings[it->a_id()].push_back(((it->a_rc() ? it->a_length() - it->a_end() : it->a_begin()) + 1) << 1 | 0);
            ref_mappings[it->a_id()].push_back(((it->a_rc() ? it->a_length() - it->a_begin() : it->a_end()) - 1) << 1 | 1);
        }
        overlaps.clear();
        hill_regions.resize(std::max(ref_mappings.size(), mappings.size()));

        std::vector<std::future<void>> thread_futures;
        for (uint32_t i = 0; i < ref_mappings.size(); ++i) {
            thread_futures.emplace_back(thread_pool->submit_task(findHillsWithReference,
                std::ref(hill_regions[i]), std::ref(ref_mappings[i]), i, false));
        }
        for (uint32_t i = 0; i < ref_mappings.size(); ++i) {
            thread_futures[i].wait();
            std::vector<uint32_t>().swap(ref_mappings[i]);
        }
        std::vector<std::vector<uint32_t>>().swap(ref_mappings);

        oreader.reset();
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

            if (hill_regions[it->a_id()].size() != 0) {
                uint32_t begin = it->a_rc() ? it->a_length() - it->a_end() : it->a_begin();
                uint32_t end = it->a_rc() ? it->a_length() - it->a_begin() : it->a_end();

                for (const auto& h: hill_regions[it->a_id()]) {
                    if (begin < h.second && h.first < end) {
                        // if (h.second < it->a_length() / 2) {
                        if (h.first < 0.10 * it->a_length()) {
                            if (end < h.second) {
                                is_valid_overlap[it->id()] = false;
                            }
                        // } else {
                        } else if (h.second > 0.9 * it->a_length()) {
                            if (begin > h.first) {
                                is_valid_overlap[it->id()] = false;
                            }
                        }
                    }
                }
            }

            if (hill_regions[it->b_id()].size() != 0) {
                uint32_t begin = it->b_rc() ? it->b_length() - it->b_end() : it->b_begin();
                uint32_t end = it->b_rc() ? it->b_length() - it->b_begin() : it->b_end();

                for (const auto& h: hill_regions[it->b_id()]) {
                    if (begin < h.second && h.first < end) {
                        // if (h.second < it->b_length() / 2) {
                        if (h.first < 0.10 * it->b_length()) {
                            if (end < h.second) {
                                is_valid_overlap[it->id()] = false;
                            }
                        // } else {
                        } else if (h.second > 0.9 * it->b_length()) {
                            if (begin > h.first) {
                                is_valid_overlap[it->id()] = false;
                            }
                        }
                    }
                }
            }

            if (prefilter) {
                bool is_valid = it->update(
                    longest_region[it->a_id()].first,
                    longest_region[it->a_id()].second,
                    longest_region[it->b_id()].first,
                    longest_region[it->b_id()].second
                );

                if (!is_valid) {
                    is_valid_overlap[it->id()] = false;
                    continue;
                }
            }

            it->set_type(classifyOverlap(it));
            switch (it->type()) {
                case 0:
                    is_valid_overlap[it->id()] = false;
                    break;
                case 1:
                    is_valid_read[it->a_id()] = false;
                    is_valid_overlap[it->id()] = false;
                    break;
                case 2:
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

    auto rreader = bioparser::createReader<Read, bioparser::FastqReader>(reads_path);

    // uint32_t trimmed_read_id = 1;
    // std::ofstream trimmed_reads_file;
    /*if (prefilter) {
        trimmed_reads_file.open("trimmed_reads.fastq");
    }*/

    while (true) {
        uint64_t current_read_id = reads.size();
        auto status = rreader->read_objects(reads, kChunkSize);

        for (uint64_t i = current_read_id; i < reads.size(); ++i) {
            auto& it = reads[i];
            // print valid read regions for later use!
            /*if (prefilter && !valid_regions[it->id()].empty()) {
                for (const auto& r: valid_regions[it->id()]) {
                    trimmed_reads_file << "@" << trimmed_read_id++ << std::endl;
                    trimmed_reads_file << it->sequence().substr(r.first, r.second - r.first) << std::endl;
                    trimmed_reads_file << "+" << std::endl;
                    trimmed_reads_file << it->quality().substr(r.first, r.second - r.first) << std::endl;
                }
            }*/

            if (it->id() >= is_valid_read.size()) {
                is_valid_read.resize(it->id() + 1, false);
                continue;
            }
            if (!is_valid_read[it->id()]) {
                continue;
            }
            if (prefilter) {
                it->trim_sequence(longest_region[it->id()].first, longest_region[it->id()].second);
            }
        }

        shrinkVector(reads, current_read_id, is_valid_read);

        if (status == false) {
            break;
        }
    }

    /*if (prefilter) {
        trimmed_reads_file.close();
    }*/

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
