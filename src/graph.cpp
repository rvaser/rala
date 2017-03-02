/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <set>
#include <list>
#include <stack>
#include <deque>
#include <tuple>
#include <algorithm>

#include <iostream>
#include <fstream>

#include "read.hpp"
#include "overlap.hpp"
#include "graph.hpp"
#include "utils.hpp"

namespace rala {

constexpr double kChimericRatio = 1.75;
constexpr double kMinMatchingBasesRatio = 2.5;
constexpr uint32_t kMinCoverage = 3;
constexpr double kMaxOverhangToOverlapRatio = 0.875; //0.8;
constexpr double kTransitiveEdgeEps = 0.12;
constexpr uint32_t kMaxBubbleLength = 5000000;
constexpr uint32_t kMinUnitigSize = 6;

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

uint32_t classifyOverlap(const std::shared_ptr<Overlap>& overlap, double& tt) {

    uint32_t left_overhang = std::min(overlap->a_begin(), overlap->b_begin());
    uint32_t right_overhang = std::min(overlap->a_length() - overlap->a_end(),
        overlap->b_length() - overlap->b_end());

    uint32_t a_len = overlap->a_end() - overlap->a_begin();
    uint32_t b_len = overlap->b_end() - overlap->b_begin();

    tt += a_len / (double) (a_len + left_overhang + right_overhang);
    tt += b_len / (double) (b_len + left_overhang + right_overhang);

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

// call after hit sorting!!
void printCoverageGraph(const std::vector<uint32_t>& hits, const std::string& path) {

    std::vector<uint32_t> coverage_graph((hits.back() >> 1) + 1, 0);
    int32_t coverage = 0;
    uint32_t last = 0, total = 0;
    for (const auto& hit: hits) {
        if (coverage > 0) {
            for (uint32_t i = last; i < (hit >> 1); ++i) {
                coverage_graph[i] += coverage;
                total += coverage;
            }
        }
        last = hit >> 1;
        if (hit & 1) --coverage;
        else ++coverage;
    }
    std::ofstream out(path);
    for (uint32_t i = 0; i < coverage_graph.size(); ++i) {
        out << i << " " << coverage_graph[i] << " " << 0 << std::endl;
    }
    out.close();
}

bool slopeRead(std::vector<std::pair<uint32_t, uint32_t>>& slope_regions, const std::vector<uint32_t>& hits,
    const std::string& path, bool print = false) {

    bool is_slope_read = false;
    std::pair<uint32_t, uint32_t> slope_region(-1, 0);
    // slope_region.first = -1;
    // slope_region.second = 0;

    auto deque_add = [](std::deque<std::pair<int32_t, int32_t>>& window, int32_t k, int32_t value, int32_t position) -> void {
        while (!window.empty() && window.back().second <= value) {
            window.pop_back();
        }
        window.emplace_back(position, value);
        return;
    };

    auto deque_update = [](std::deque<std::pair<int32_t, int32_t>>& window, int32_t k, int32_t position) -> void {
        while (!window.empty() && window.front().first <= position - k) {
            window.pop_front();
        }
        return;
    };

    std::vector<uint32_t> coverage_graph((hits.back() >> 1) + 1, 0);
    int32_t coverage = 0;
    uint32_t last = 0;
    for (const auto& hit: hits) {
        if (coverage > 0) {
            for (uint32_t i = last; i < (hit >> 1); ++i) {
                coverage_graph[i] += coverage;
            }
        }
        last = hit >> 1;
        if (hit & 1) --coverage;
        else ++coverage;
    }

    std::deque<std::pair<int32_t, int32_t>> left_window, right_window;
    std::vector<int32_t> slopes;

    int32_t k = std::max(500U, uint32_t(0.03 * ((hits.back()>>1) - (hits.front()>>1))));
    int32_t length = coverage_graph.size();

    for (int32_t i = -1 * k + 2; i < length - 1; ++i) {
        if (i < length - k) {
            deque_add(right_window, k, coverage_graph[i + k], i + k);
        }
        deque_update(right_window, k, i + k);

        if (i > 0) {
            deque_add(left_window, k, coverage_graph[i - 1], i - 1);
            deque_update(left_window, k, i - 1);

            int32_t current = coverage_graph[i] * kChimericRatio;
            if (left_window.front().second > current) {
                slopes.push_back(i << 1 | 0);
            }
            if (right_window.front().second > current) {
                slopes.push_back(i << 1 | 1);
            }
        }
    }

    if (slopes.size() > 0) {
        std::sort(slopes.begin(), slopes.end());
        for (uint32_t i = 0; i < slopes.size() - 1; ++i) {
            if (!(slopes[i] & 1) && (slopes[i + 1] & 1)) {
                is_slope_read = true;
                slope_region.first = std::min(slope_region.first, (uint32_t) (slopes[i] >> 1));
                slope_region.second = std::max(slope_region.second, (uint32_t) (slopes[i + 1] >> 1));
            }
        }
    }

    slope_regions.push_back(slope_region);

    if (print && is_slope_read) {
        std::ofstream out(path);
        std::vector<int32_t> slopes(coverage_graph.size(), 0);
        slopes[slope_region.first] = 1;
        slopes[slope_region.second] = 2;
        for (uint32_t i = 0; i < coverage_graph.size(); ++i) {
            out << i << " " << coverage_graph[i] << " " <<  slopes[i] << std::endl;
        }
        out.close();
        fprintf(stderr, "%s (k = %d): %d - %d\n", path.c_str(), k, slope_region.first, slope_region.second);
    }

    return is_slope_read;
}

void findSlopes(std::vector<std::pair<uint32_t, uint32_t>>& hills, const std::vector<uint32_t>& hits,
    const std::string& path, bool print = false) {

    std::vector<uint32_t> coverage_graph((hits.back() >> 1) + 1, 0);
    int32_t coverage = 0;
    uint32_t last = 0;
    for (const auto& hit: hits) {
        if (coverage > 0) {
            for (uint32_t i = last; i < (hit >> 1); ++i) {
                coverage_graph[i] += coverage;
            }
        }
        last = hit >> 1;
        if (hit & 1) --coverage;
        else ++coverage;
    }

    auto window_add = [](std::deque<std::pair<int32_t, int32_t>>& window, int32_t k, int32_t value, int32_t position) -> void {
        while (!window.empty() && window.back().second <= value) {
            window.pop_back();
        }
        window.emplace_back(position, value);
        return;
    };

    auto window_update = [](std::deque<std::pair<int32_t, int32_t>>& window, int32_t k, int32_t position) -> void {
        while (!window.empty() && window.front().first <= position - k) {
            window.pop_front();
        }
        return;
    };

    std::deque<std::pair<int32_t, int32_t>> left_window, right_window;
    std::vector<int32_t> slopes;

    int32_t k = std::max(500U, uint32_t(0.03 * ((hits.back()>>1) - (hits.front()>>1))));
    int32_t read_length = coverage_graph.size();

    for (int32_t i = -1 * k + 2; i < read_length; ++i) {
        if (i < read_length - k) {
            window_add(right_window, k, coverage_graph[i + k], i + k);
        }
        window_update(right_window, k, i + k);

        if (i > 0) {
            window_add(left_window, k, coverage_graph[i - 1], i - 1);
            window_update(left_window, k, i - 1);

            int32_t current = coverage_graph[i] * 1.3;
            if (left_window.front().second > current) {
                slopes.push_back(i << 1 | 0);
            }
            if (!right_window.empty() && right_window.front().second > current) {
                slopes.push_back(i << 1 | 1);
            }
        }
    }

    int32_t max_width = ((hits.back()>>1) - (hits.front()>>1)) * 0.7;

    if (slopes.size() > 2) {

        fprintf(stderr, "%s\n", path.c_str());
        uint32_t ldownslope = 0, fdownslope = 0;
        bool found_fds = false;

        uint32_t fupslope = 0, lupslope = 0;
        bool found_fus = false;

        std::vector<std::pair<uint32_t, uint32_t>> slope_regions;

        uint32_t slope_width = 2 * k;
        for (uint32_t s = 1; s < slopes.size(); ++s) {
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

        for (const auto& it: slope_regions) {
            fprintf(stderr, "%d (%d %d) \n", it.first & 1, it.first >> 1, it.second);
        }

        // trim regions
        auto trim_left = [&](uint32_t begin, uint32_t end) -> uint32_t {
            uint32_t value = coverage_graph[begin];
            while (begin <= end && coverage_graph[begin] == value) {
                ++begin;
            }
            return begin - 1;
        };

        auto trim_right = [&](uint32_t begin, uint32_t end) -> uint32_t {
            uint32_t value = coverage_graph[end];
            while (end >= begin && coverage_graph[end] == value) {
                --end;
            }
            return end + 1;
        };

        for (auto& region: slope_regions) {
            if (region.first & 1) { // up slope
                region.first = (trim_left(region.first >> 1, region.second) << 1) | 1;
                region.second = trim_right(region.first >> 1, region.second);
            } else { // downslope
                region.second = trim_right(region.first >> 1, region.second);
                region.first = (trim_left(region.first >> 1, region.second) << 1) | 0;
            }
        }
        fprintf(stderr, "TRIMMED\n");
        for (const auto& it: slope_regions) {
            fprintf(stderr, "%d (%d %d) \n", it.first & 1, it.first >> 1, it.second);
        }

        // chimeric check
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
        fprintf(stderr, "Chimeric check\n");
        for (const auto& it: slope_regions) {
            fprintf(stderr, "%d (%d %d) \n", it.first & 1, it.first >> 1, it.second);
        }
        fprintf(stderr, "\n");

        std::vector<uint32_t> slope_graph(coverage_graph.size(), 0);
        slopes.clear();
        for (const auto& region: slope_regions) {
            if ((region.first >> 1) < (hits.front() >> 1)) {
                // slopes.push_back((((region.first >> 1) + region.second) / 2 + 1) << 1 | (region.first & 1));
                slopes.push_back((region.second + 1) << 1 | 1);
            } else if (region.second >= (hits.back() >> 1)) {
                // slopes.push_back((((region.first >> 1) + region.second) / 2 - 1) << 1 | (region.first & 1));
                slopes.push_back(((region.first >> 1) - 1) << 1 | 0);
            } else if (region.first & 1) { // upslope
                slopes.push_back(region.first);
            } else { // downslope
                slopes.push_back(region.second << 1 | 0);
            }
            // slope_graph[slopes.back() >> 1] = slopes.back() & 1 ? 2 : 1;
        }

        // find hills
        bool has_hill = false;
        for (uint32_t s = 0; s < slopes.size() - 1; ++s) {
            if ((slopes[s] & 1) && !(slopes[s + 1] & 1)) {
                if ((slopes[s + 1] >> 1) - (slopes[s] >> 1) < max_width) {
                    slope_graph[(slopes[s] >> 1)] = 2;
                    slope_graph[(slopes[s + 1] >> 1)] = 1;
                    has_hill = true;
                    hills.emplace_back(slopes[s] >> 1, slopes[s + 1] >> 1);
                }
            }
        }

        if (has_hill && print) {
            std::ofstream out(path);
            for (uint32_t i = 0; i < coverage_graph.size(); ++i) {
                out << i << " " << coverage_graph[i] << " " << slope_graph[i] << std::endl;
            }
            out.close();
        }
    }
}

// call before hit sorting!!
void printPileGraph(const std::vector<uint32_t>& hits, const std::string& path) {

    std::ofstream out(path);
    for (uint32_t i = 0; i < hits.size(); i+= 2) {
        out << (hits[i] >> 1) << " " << (hits[i+1] >> 1) << std::endl;
    }
    out.close();
};

void prefilterData(std::vector<std::pair<uint32_t, uint32_t>>& regions,
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& slope_regions,
    std::vector<bool>& is_valid_read,
    std::vector<bool>& is_valid_overlap,
    const std::string& reads_path,
    const std::string& overlaps_path,
    uint32_t overlap_type) {

    fprintf(stderr, "  Prefiltering data {\n");

    std::vector<std::shared_ptr<Overlap>> overlaps;
    auto oreader = overlap_type == 0 ?
        bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
        bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

    std::vector<std::vector<uint32_t>> hits(is_valid_read.size());

    while (true) {
        auto status = oreader->read_objects(overlaps, kChunkSize);

        for (const auto& it: overlaps) {
            if (!is_valid_overlap[it->id()]) {
                continue;
            }

            uint32_t begin = (it->a_rc() ? it->a_length() - it->a_end() : it->a_begin()) + 1;
            uint32_t end = (it->a_rc() ? it->a_length() - it->a_begin() : it->a_end()) - 1;
            hits[it->a_id()].push_back(begin << 1 | 0);
            hits[it->a_id()].push_back(end << 1 | 1);

            begin = (it->b_rc() ? it->b_length() - it->b_end() : it->b_begin()) + 1;
            end = (it->b_rc() ? it->b_length() - it->b_begin() : it->b_end()) - 1;
            hits[it->b_id()].push_back(begin << 1 | 0);
            hits[it->b_id()].push_back(end << 1 | 1);
        }

        overlaps.clear();

        if (status == false) {
            break;
        }
    }
    oreader.reset();

    uint32_t num_slope_reads = 0;
    uint32_t num_trimmed_reads = 0;

    for (uint32_t i = 0; i < hits.size(); ++i) {
        if (hits[i].empty()) {
            is_valid_read[i] = false;
            continue;
        }

        // if (i == 24337) printPileGraph(hits[i], "graphs/" + std::to_string(i));
        std::sort(hits[i].begin(), hits[i].end());

        findSlopes(slope_regions[i], hits[i], "graphs/l" + std::to_string(i));
        //if (slopeRead(slope_regions[i], hits[i], "graphs/c" + std::to_string(i))) {
        //    ++num_slope_reads;
        //}

        int32_t start = hits[i].front() >> 1;
        int32_t end = hits[i].back() >> 1;

        if (start > end) {
            continue;
        }

        int32_t coverage = 0, min_coverage = kMinCoverage;
        int32_t begin = 0, max_begin = 0, max_end = 0;
        for (const auto& hit: hits[i]) {
            int32_t old_coverage = coverage;
            if (hit & 1) {
                --coverage;
            } else {
                ++coverage;
            }
            int32_t pos = hit >> 1;
            if (old_coverage < min_coverage && coverage >= min_coverage) {
                begin = pos;
            } else if (old_coverage >= min_coverage && coverage < min_coverage) {
                int32_t length = pos - begin;
                if (length > max_end - max_begin) {
                    max_begin = begin;
                    max_end = pos;
                }
            }
        }

        if (max_end - max_begin > 0) {
            regions[i] = std::make_pair(max_begin, max_end);
            ++num_trimmed_reads;
        } else {
            is_valid_read[i] = false;
        }

        std::vector<uint32_t>().swap(hits[i]);
    }
    std::vector<std::vector<uint32_t>>().swap(hits);

    // save trimmed fastq reads for consensus!

    fprintf(stderr, "    number of slope reads = %u\n", num_slope_reads);
    fprintf(stderr, "    number of trimmed reads = %u\n", num_trimmed_reads);
    fprintf(stderr, "  }\n");
}

void preprocessData(std::vector<std::shared_ptr<Read>>& reads, std::vector<std::shared_ptr<Overlap>>& overlaps,
    const std::string& reads_path, const std::string& overlaps_path, uint32_t overlap_type, bool prefilter) {

    fprintf(stderr, "Preprocessing data {\n");

    auto remove_multiple_overlaps = [](const std::vector<std::shared_ptr<Overlap>>& overlaps,
        std::vector<bool>& is_valid_overlap) -> uint32_t {

        uint32_t num_multiple_matches = 0;
        for (uint32_t i = 0; i < overlaps.size(); ++i) {
            if (is_valid_overlap[overlaps[i]->id()] == false) {
                continue;
            }
            for (uint32_t j = i + 1; j < overlaps.size(); ++j) {
                if (is_valid_overlap[overlaps[j]->id()] == false ||
                    overlaps[i]->b_id() != overlaps[j]->b_id()) {
                    continue;
                }
                ++num_multiple_matches;
                if (overlaps[i]->length() > overlaps[j]->length()) {
                    is_valid_overlap[overlaps[i]->id()] = false;
                    break;
                } else {
                    is_valid_overlap[overlaps[j]->id()] = false;
                }
            }
        }

        return num_multiple_matches;
    };

    std::vector<bool> is_valid_overlap;
    std::vector<bool> is_valid_read;
    uint32_t max_read_id = 0;

    std::vector<std::shared_ptr<Overlap>> current_overlaps;
    uint32_t num_multiple_matches = 0;
    uint32_t num_self_matches = 0;

    auto oreader = overlap_type == 0 ?
        bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
        bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

    while (true) {

        auto status = oreader->read_objects(overlaps, kChunkSize);
        is_valid_overlap.resize(is_valid_overlap.size() + overlaps.size(), true);

        for (const auto& it: overlaps) {

            max_read_id = std::max(max_read_id, std::max(it->a_id(), it->b_id()));

            // self match check
            if (it->a_id() == it->b_id()) {
                is_valid_overlap[it->id()] = false;
                ++num_self_matches;
                continue;
            }

            // multiple overlaps check
            if (current_overlaps.size() != 0 && current_overlaps.front()->a_id() != it->a_id()) {
                num_multiple_matches += remove_multiple_overlaps(current_overlaps, is_valid_overlap);
                current_overlaps.clear();
            }
            current_overlaps.push_back(it);
        }

        overlaps.clear();

        if (status == false) {
            num_multiple_matches += remove_multiple_overlaps(current_overlaps, is_valid_overlap);
            current_overlaps.clear();
            break;
        }
    }
    oreader.reset();

    is_valid_read.resize(max_read_id + 1, true);

    std::vector<std::pair<uint32_t, uint32_t>> regions(is_valid_read.size());
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> slope_regions(is_valid_read.size());
    if (prefilter == true) {
        prefilterData(regions, slope_regions, is_valid_read, is_valid_overlap,
            reads_path, overlaps_path, overlap_type);
    }

    // reading valid overlaps into memory
    oreader = overlap_type == 0 ?
        bioparser::createReader<Overlap, bioparser::MhapReader>(overlaps_path) :
        bioparser::createReader<Overlap, bioparser::PafReader>(overlaps_path);

    std::vector<uint32_t> num_overlaps_per_read(is_valid_read.size(), 0);
    std::vector<uint32_t> num_removed_hillaps(is_valid_read.size(), 0);

    uint32_t num_slope_overlaps = 0;
    std::set<uint32_t> slope_read_overlaps;
    uint64_t tot = 0;
    double tt = 0;
    while (true) {
        uint64_t current_overlap_id = overlaps.size();
        // uint64_t i = start;
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
            ++num_overlaps_per_read[it->a_id()];
            ++num_overlaps_per_read[it->b_id()];

            if (prefilter) {
                /*if (slope_regions[it->a_id()].front().second != 0) {
                    slope_read_overlaps.insert(it->id());
                    uint32_t begin = (it->a_rc() ? it->a_length() - it->a_end() : it->a_begin());
                    uint32_t end = (it->a_rc() ? it->a_length() - it->a_begin() : it->a_end());
                    if (begin > slope_regions[it->a_id()].front().first || end < slope_regions[it->a_id()].front().second) {
                        is_valid_overlap[it->id()] = false;
                        ++num_slope_overlaps;
                        continue;
                    }
                }
                if (slope_regions[it->b_id()].front().second != 0) {
                    slope_read_overlaps.insert(it->id());
                    uint32_t begin = (it->b_rc() ? it->b_length() - it->b_end() : it->b_begin());
                    uint32_t end = (it->b_rc() ? it->b_length() - it->b_begin() : it->b_end());
                    if (begin > slope_regions[it->b_id()].front().first || end < slope_regions[it->b_id()].front().second) {
                        is_valid_overlap[it->id()] = false;
                        ++num_slope_overlaps;
                        continue;
                    }
                }*/
                uint32_t contained_hills = 0, overlaping_hills = 0;
                if (slope_regions[it->a_id()].size() != 0) {
                    slope_read_overlaps.insert(it->id());
                    uint32_t begin = it->a_rc() ? it->a_length() - it->a_end() : it->a_begin();
                    uint32_t end = it->a_rc() ? it->a_length() - it->a_begin() : it->a_end();

                    for (const auto& sreg: slope_regions[it->a_id()]) {
                        if (begin < sreg.first && end > sreg.second) {
                            ++contained_hills;
                        } else if (!(begin > sreg.second || end < sreg.first)) {
                            ++overlaping_hills;
                        }
                    }
                }

                if (slope_regions[it->b_id()].size() != 0) {
                    slope_read_overlaps.insert(it->id());
                    uint32_t begin = it->b_rc() ? it->b_length() - it->b_end() : it->b_begin();
                    uint32_t end = it->b_rc() ? it->b_length() - it->b_begin() : it->b_end();

                    for (const auto& sreg: slope_regions[it->b_id()]) {
                        if (begin < sreg.first && end > sreg.second) {
                            ++contained_hills;
                        } else if (!(begin > sreg.second || end < sreg.first)) {
                            ++overlaping_hills;
                        }
                    }
                }

                if (overlaping_hills > 0 && contained_hills == 0) {
                    ++num_removed_hillaps[it->a_id()];
                    ++num_removed_hillaps[it->b_id()];
                    ++num_slope_overlaps;
                    is_valid_overlap[it->id()] = false;
                    // continue;
                }

                bool is_valid = it->update(
                    regions[it->a_id()].first,
                    regions[it->a_id()].second,
                    regions[it->b_id()].first,
                    regions[it->b_id()].second
                );

                if (!is_valid) {
                    is_valid_overlap[it->id()] = false;
                    continue;
                }
            }

            it->set_type(classifyOverlap(it, tt));
            tot += 2;
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

    uint32_t zero_ovl_reads = 0;
    for (uint32_t i = 0; i < is_valid_read.size(); ++i) {
        if (num_overlaps_per_read[i] != 0 && num_overlaps_per_read[i] == num_removed_hillaps[i]) {
            // fprintf(stderr, "0ovl read = %d\n", i);
            ++zero_ovl_reads;
        }
    }
    fprintf(stderr, "Zero ovl reads = %d\n", zero_ovl_reads);
    fprintf(stderr, "Avg len/ovh ratio = %lf, %lf, %lu\n", tt / (double) tot, tt, tot);

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

    while (true) {
        uint64_t current_read_id = reads.size();
        // uint64_t r = start;
        auto status = rreader->read_objects(reads, kChunkSize);

        for (uint64_t i = current_read_id; i < reads.size(); ++i) {
            auto& it = reads[i];
            if (it->id() >= is_valid_read.size()) {
                is_valid_read.resize(it->id() + 1, false);
                continue;
            }
            if (!is_valid_read[it->id()]) {
                continue;
            }
            if (prefilter) {
                it->trim_sequence(regions[it->id()].first, regions[it->id()].second);
            }
        }

        shrinkVector(reads, current_read_id, is_valid_read);

        if (status == false) {
            break;
        }
    }
    rreader.reset();

    uint64_t num_valid_reads = 0, num_valid_overlaps = 0;
    for (const auto& it: is_valid_read) if (it == true) ++num_valid_reads;
    for (const auto& it: is_valid_overlap) if (it == true) ++num_valid_overlaps;

    fprintf(stderr, "  number of self matches = %u\n", num_self_matches);
    fprintf(stderr, "  number of multiple matches = %u\n", num_multiple_matches);
    fprintf(stderr, "  number of slope overlaps = %u (out of %lu)\n", num_slope_overlaps, slope_read_overlaps.size());
    fprintf(stderr, "  number of valid overlaps = %lu (out of %lu)\n", num_valid_overlaps, is_valid_overlap.size());
    fprintf(stderr, "  number of valid reads = %lu (out of %lu)\n", num_valid_reads, is_valid_read.size());
    fprintf(stderr, "}\n\n");
}

class Graph::Node {
    public:
        // Node encapsulating read
        Node(uint32_t _id, const std::shared_ptr<Read>& read) :
                id(_id), read_id(read->id()), pair(), sequence(id % 2 == 0 ? read->sequence() : read->rc()),
                prefix_edges(), suffix_edges(), read_ids(1, read->id()), unitig_size(1), mark(false) {

            auto is_unitig = read->name().find("Utg=");
            if (is_unitig != std::string::npos) {
                unitig_size = atoi(read->name().c_str() + is_unitig + 4);
                // fprintf(stderr, "Unitig size = %u\n", unitig_size);
            }
        }
        // Unitig
        Node(uint32_t _id, Node* begin_node, Node* end_node, std::unordered_set<uint32_t>& marked_edges);
        // Circular unitig
        Node(uint32_t _id, Node* begin_node, std::unordered_set<uint32_t>& marked_edges);
        Node(const Node&) = delete;
        const Node& operator=(const Node&) = delete;

        ~Node() {}

        uint32_t length() const {
            return sequence.size();
        }

        uint32_t in_degree() const {
            return prefix_edges.size();
        }

        uint32_t out_degree() const {
            return suffix_edges.size();
        }

        bool is_junction() const {
            return (out_degree() > 1 || in_degree() > 1);
        }

        bool is_tip() const {
            return (out_degree() > 0 && in_degree() == 0 && unitig_size < kMinUnitigSize);
        }

        uint32_t id;
        uint32_t read_id;
        Node* pair;
        std::string sequence;
        std::list<Edge*> prefix_edges;
        std::list<Edge*> suffix_edges;
        std::vector<uint32_t> read_ids;
        uint32_t unitig_size;
        bool mark;
};

class Graph::Edge {
    public:
        Edge(uint32_t _id, const std::shared_ptr<Overlap>& overlap, Node* _begin_node,
            Node* _end_node, uint32_t type) :
                id(_id), pair(), begin_node(_begin_node), end_node(_end_node), length(),
                quality(overlap->quality()), mark(false) {

            uint32_t length_a = id % 2 == 0 ? overlap->a_begin() : overlap->a_length() - overlap->a_end();
            uint32_t length_b = id % 2 == 0 ? overlap->b_begin() : overlap->b_length() - overlap->b_end();

            if (type == 0) { // a to b overlap
                length = length_a - length_b;
            } else { // b to a overlap
                length = length_b - length_a;
            }
        }
        Edge(const Edge&) = delete;
        const Edge& operator=(const Edge&) = delete;

        ~Edge() {}

        std::string label() const {
            return begin_node->sequence.substr(0, length);
        }

        uint32_t matching_bases() const {
            return (quality * (begin_node->length() - length));
        }

        uint32_t id;
        Edge* pair;
        Node* begin_node;
        Node* end_node;
        uint32_t length;
        double quality;
        bool mark;
};

Graph::Node::Node(uint32_t _id, Node* begin_node, Node* end_node, std::unordered_set<uint32_t>& marked_edges) :
        id(_id), read_id(), pair(), sequence(), prefix_edges(), suffix_edges(),
        read_ids(), unitig_size(), mark(false) {

    if (!begin_node->prefix_edges.empty()) {
        begin_node->prefix_edges.front()->end_node = this;
        prefix_edges.push_back(begin_node->prefix_edges.front());
    }

    uint32_t length = 0;
    Node* curr_node = begin_node;
    while (curr_node->id != end_node->id) {
        auto* edge = curr_node->suffix_edges.front();
        edge->mark = true;
        marked_edges.insert(edge->id);

        read_ids.reserve(read_ids.size() + curr_node->read_ids.size());
        read_ids.insert(read_ids.end(), curr_node->read_ids.begin(), curr_node->read_ids.end());

        unitig_size += curr_node->unitig_size;
        length += edge->length;
        sequence += edge->label();

        curr_node->prefix_edges.clear();
        curr_node->suffix_edges.clear();
        curr_node->mark = true;

        curr_node = edge->end_node;
    }

    read_ids.reserve(read_ids.size() + end_node->read_ids.size());
    read_ids.insert(read_ids.end(), end_node->read_ids.begin(), end_node->read_ids.end());

    unitig_size += end_node->unitig_size;
    sequence += end_node->sequence;

    if (!end_node->suffix_edges.empty()) {
        end_node->suffix_edges.front()->begin_node = this;
        end_node->suffix_edges.front()->length += length;
        suffix_edges.push_back(end_node->suffix_edges.front());
    }

    end_node->prefix_edges.clear();
    end_node->suffix_edges.clear();
    end_node->mark = true;
}

Graph::Node::Node(uint32_t _id, Node* begin_node, std::unordered_set<uint32_t>& marked_edges) :
        id(_id), read_id(), pair(), sequence(), prefix_edges(), suffix_edges(),
        read_ids(), unitig_size(), mark(false) {
    // fprintf(stderr, "!!! CIRCULAR UNITIG ALERT !!!\n");

    uint32_t length = 0;
    Node* curr_node = begin_node;
    while (true) {
        // fprintf(stderr, "Curr node = %d\n", curr_node->id);
        auto* edge = curr_node->suffix_edges.front();
        edge->mark = true;
        marked_edges.insert(edge->id);

        read_ids.reserve(read_ids.size() + curr_node->read_ids.size());
        read_ids.insert(read_ids.end(), curr_node->read_ids.begin(), curr_node->read_ids.end());

        unitig_size += curr_node->unitig_size;
        length += edge->length;
        sequence += edge->label();

        curr_node->prefix_edges.clear();
        curr_node->suffix_edges.clear();
        curr_node->mark = true;

        curr_node = edge->end_node;
        if (curr_node->id == begin_node->id) {
            break;
        }
    }
}

std::unique_ptr<Graph> createGraph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps) {
    return std::unique_ptr<Graph>(new Graph(reads, overlaps));
}

Graph::Graph(const std::vector<std::shared_ptr<Read>>& reads,
    const std::vector<std::shared_ptr<Overlap>>& overlaps)
        : nodes_(), edges_(), marked_edges_() {

    fprintf(stderr, "Assembly graph {\n");

    uint64_t max_read_id = 0;
    for (const auto& it: reads) {
        max_read_id = std::max(max_read_id, it->id());
    }

    // create assembly graph
    std::vector<int32_t> read_id_to_node_id(max_read_id + 1, -1);
    uint32_t node_id = 0;
    for (const auto& read: reads) {
        read_id_to_node_id[read->id()] = node_id;

        Node* node = new Node(node_id++, read); // normal read
        Node* _node = new Node(node_id++, read); // reverse complement

        node->pair = _node;
        _node->pair = node;

        nodes_.push_back(std::unique_ptr<Node>(node));
        nodes_.push_back(std::unique_ptr<Node>(_node));
    }

    uint32_t edge_id = 0;
    for (const auto& overlap: overlaps) {

        auto a = nodes_[read_id_to_node_id[overlap->a_id()] + (overlap->a_rc() == 0 ? 0 : 1)].get();
        auto _a = a->pair;

        auto b = nodes_[read_id_to_node_id[overlap->b_id()] + (overlap->b_rc() == 0 ? 0 : 1)].get();
        auto _b = b->pair;

        if (overlap->type() == 3) { // a to b overlap
            Edge* edge = new Edge(edge_id++, overlap, a, b, 0);
            Edge* _edge = new Edge(edge_id++, overlap, _b, _a, 1);

            edge->pair = _edge;
            _edge->pair = edge;

            edges_.push_back(std::unique_ptr<Edge>(edge));
            edges_.push_back(std::unique_ptr<Edge>(_edge));

            a->suffix_edges.push_back(edge);
            _a->prefix_edges.push_back(_edge);
            b->prefix_edges.push_back(edge);
            _b->suffix_edges.push_back(_edge);

        } else if (overlap->type() == 4) { // b to a overlap
            Edge* edge = new Edge(edge_id++, overlap, b, a, 1);
            Edge* _edge = new Edge(edge_id++, overlap, _a, _b, 0);

            edge->pair = _edge;
            _edge->pair = edge;

            edges_.push_back(std::unique_ptr<Edge>(edge));
            edges_.push_back(std::unique_ptr<Edge>(_edge));

            b->suffix_edges.push_back(edge);
            _b->prefix_edges.push_back(_edge);
            a->prefix_edges.push_back(edge);
            _a->suffix_edges.push_back(_edge);
        }
    }

    fprintf(stderr, "  Construction {\n");
    fprintf(stderr, "    number of graph nodes = %zu\n", nodes_.size());
    fprintf(stderr, "    number of graph edges = %zu\n", edges_.size());
    fprintf(stderr, "  }\n");
}

Graph::~Graph() {
    fprintf(stderr, "}\n");
}

void Graph::remove_isolated_nodes() {

    for (auto& node: nodes_) {
        if (node == nullptr) {
            continue;
        }
        if ((node->in_degree() == 0 && node->out_degree() == 0 && node->unitig_size < kMinUnitigSize) || (node->mark == true)) {
            // fprintf(stderr, "Removing isolated node: %d\n", node->id);
            node.reset();
        }
    }
}

void Graph::remove_transitive_edges() {

    fprintf(stderr, "  Transitive edge removal {\n");

    uint32_t num_transitive_edges = 0;
    std::vector<Edge*> candidate_edge(nodes_.size(), nullptr);

    for (const auto& node_x: nodes_) {
        if (node_x == nullptr) continue;

        for (const auto& edge: node_x->suffix_edges) {
            candidate_edge[edge->end_node->id] = edge;
        }

        for (const auto& edge_xy: node_x->suffix_edges) {
            for (const auto& edge_yz: nodes_[edge_xy->end_node->id]->suffix_edges) {
                uint32_t z = edge_yz->end_node->id;
                if (candidate_edge[z] != nullptr && candidate_edge[z]->mark == false) {
                    if (isSimilar(edge_xy->length + edge_yz->length, candidate_edge[z]->length, kTransitiveEdgeEps)) {
                        candidate_edge[z]->mark = true;
                        candidate_edge[z]->pair->mark = true;
                        marked_edges_.insert(candidate_edge[z]->id);
                        marked_edges_.insert(candidate_edge[z]->pair->id);
                        ++num_transitive_edges;
                    }
                }
            }
        }

        for (const auto& edge: node_x->suffix_edges) {
            candidate_edge[edge->end_node->id] = nullptr;
        }
    }
    remove_marked_edges();

    fprintf(stderr, "    removed %u edges\n", num_transitive_edges);
    fprintf(stderr, "  }\n");
}

void Graph::remove_long_edges() {

    fprintf(stderr, "  Long edge removal {\n");
    uint32_t num_long_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || node->suffix_edges.size() < 2) continue;

        for (const auto& edge1: node->suffix_edges) {
            for (const auto& edge2: node->suffix_edges) {
                if (edge1->id == edge2->id || edge1->mark == true || edge2->mark == true) continue;
                if (edge1->matching_bases() > kMinMatchingBasesRatio * edge2->matching_bases()) {
                    edge2->mark = true;
                    edge2->pair->mark = true;
                    marked_edges_.insert(edge2->id);
                    marked_edges_.insert(edge2->pair->id);
                    ++num_long_edges;
                }
            }
        }
    }
    // fprintf(stderr, "Number of alignments to be done = %u\n", knots);

    fprintf(stderr, "    removed %u edges\n", num_long_edges);
    fprintf(stderr, "  }\n");

    remove_marked_edges();
}

uint32_t Graph::remove_tips() {

    fprintf(stderr, "  Tip removal {\n");

    uint32_t num_tip_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr) {
            continue;
        }
        // fprintf(stderr, "Considering node %d for tip removal\r", node->id);
        if (!node->is_tip()) {
            continue;
        }

        uint32_t num_removed_edges = 0;

        for (const auto& edge: node->suffix_edges) {
            if (edge->end_node->in_degree() > 1) {
                // fprintf(stderr, "Removing %d\n", edge->begin_node->id);
                edge->mark = true;
                edge->pair->mark = true;
                marked_edges_.insert(edge->id);
                marked_edges_.insert(edge->pair->id);
                ++num_removed_edges;
            }
        }

        if (num_removed_edges == node->suffix_edges.size()) {
            node->mark = true;
            node->pair->mark = true;
        }

        num_tip_edges += num_removed_edges;

        remove_marked_edges();
    }
    remove_isolated_nodes();

    fprintf(stderr, "    removed %u edges\n", num_tip_edges);
    fprintf(stderr, "  }\n");

    return num_tip_edges;
}

void Graph::remove_cycles() {

    fprintf(stderr, "  Cycle removal {\n");

    std::stack<uint32_t> stack;
    std::vector<int32_t> indexes(nodes_.size(), -1);
    std::vector<int32_t> low_links(nodes_.size(), -1);
    std::vector<bool> is_on_stack(nodes_.size(), false);
    int32_t index = 0;

    std::vector<std::vector<uint32_t>> cycles;

    std::function<void(uint32_t)> strong_connect = [&](uint32_t v) -> void {
        indexes[v] = index;
        low_links[v] = index;
        ++index;
        // fprintf(stderr, "Pushing %d\n", v);
        stack.push(v);
        is_on_stack[v] = true;

        for (const auto& edge: nodes_[v]->suffix_edges) {
            uint32_t w = edge->end_node->id;
            if (indexes[w] == -1) {
                strong_connect(w);
                low_links[v] = std::min(low_links[v], low_links[w]);
            } else if (is_on_stack[w]) {
                low_links[v] = std::min(low_links[v], indexes[w]);
            }
        }

        if (low_links[v] == indexes[v]) {
            // new strongly connected component
            std::vector<uint32_t> scc = { v };
            uint32_t w;
            do {
                w = stack.top();
                stack.pop();
                is_on_stack[w] = false;
                scc.push_back(w);
            } while (v != w);

            if (scc.size() > 2) {
                cycles.push_back(scc);
            }
        }
    };

    uint32_t num_cycle_edges = 0;
    do {
        cycles.clear();
        for (const auto& node: nodes_) {
            if (node == nullptr) continue;
            if (indexes[node->id] == -1) {
                strong_connect(node->id);
            }
        }

        // fprintf(stderr, "Number of cycles %zu\n", cycles.size());

        for (const auto& cycle: cycles) {

            Edge* worst_edge = nullptr;
            double min_score = 5;

            for (uint32_t i = 0; i < cycle.size() - 1; ++i) {
                const auto& node = nodes_[cycle[i]];
                for (auto& edge: node->prefix_edges) {
                    if (edge->begin_node->id == cycle[i + 1]) {
                        if (min_score > edge->quality) {
                            min_score = edge->quality;
                            worst_edge = edge;
                        }
                        break;
                    }
                }
            }

            worst_edge->mark = true;
            worst_edge->pair->mark = true;
            marked_edges_.insert(worst_edge->id);
            marked_edges_.insert(worst_edge->pair->id);
            ++num_cycle_edges;
        }

        remove_marked_edges();

    } while (cycles.size() != 0);

    fprintf(stderr, "    removed %u edges\n", num_cycle_edges);
    fprintf(stderr, "  }\n");
}

uint32_t Graph::remove_chimeras() {

    fprintf(stderr, "  Chimera removal {\n");

    std::vector<uint32_t> visited(nodes_.size(), 0);
    uint32_t visited_length = 0;
    std::vector<int32_t> predecessor(nodes_.size(), -1);
    std::deque<uint32_t> queue;

    auto extract_path = [&](std::vector<int32_t>& dst, int32_t sink) -> void {
        int32_t curr_id = sink;
        while (curr_id != -1) {
            dst.push_back(curr_id);
            curr_id = predecessor[curr_id];
        }
    };

    auto find_edge = [&](uint32_t src, uint32_t dst) -> int32_t {
        for (const auto& edge: nodes_[src]->suffix_edges) {
            if (edge->end_node->id == dst) {
                return edge->id;
            }
        }
        return -1;
    };

    auto check_path = [&](std::vector<uint32_t>& dst, const std::vector<int32_t>& path) -> void {
        // find first node with multiple in edges
        int32_t pref = -1;
        for (uint32_t i = 1; i < path.size() - 1; ++i) {
            if (nodes_[path[i]]->in_degree() > 1) {
                pref = i;
                break;
            }
        }
        // find last node with multiple out edges
        int32_t suff = -1;
        for (uint32_t i = 1; i < path.size() - 1; ++i) {
            if (nodes_[path[i]]->out_degree() > 1) {
                suff = i;
            }
        }

        if (pref == -1 && suff == -1) {
            // remove whole path
            for (uint32_t i = 0; i < path.size() - 1; ++i) {
                dst.push_back(find_edge(path[i], path[i + 1]));
            }
            return;
        }
        if (pref != -1 && nodes_[path[pref]]->out_degree() > 1) return;
        if (suff != -1 && nodes_[path[suff]]->in_degree() > 1) return;

        if (pref == -1) {
            // remove everything after last suff node
            for (uint32_t i = suff; i < path.size() - 1; ++i) {
                dst.push_back(find_edge(path[i], path[i + 1]));
            }
        } else if (suff == -1) {
            // remove everything before first pref node
            for (int32_t i = 0; i < pref; ++i) {
                dst.push_back(find_edge(path[i], path[i + 1]));
            }
        } else if (suff < pref) {
            // remove everything between last suff and first pref node
            for (int32_t i = suff; i < pref; ++i) {
                dst.push_back(find_edge(path[i], path[i + 1]));
            }
        } else if (suff >= pref && nodes_[path[0]]->in_degree() == 0) {
            // remove everything after last suff node
            for (uint32_t i = suff; i < path.size() - 1; ++i) {
                dst.push_back(find_edge(path[i], path[i + 1]));
            }
            // remove everything before first pref node
            for (int32_t i = 0; i < pref; ++i) {
                dst.push_back(find_edge(path[i], path[i + 1]));
            }
        }
    };

    uint32_t num_chimeric_edges = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || node->out_degree() < 2) {
            continue;
        }

        bool found_chim_sink = false;
        int32_t chim_sink = -1, chim_sink_pair = -1;
        int32_t source = node->id;
        std::vector<int32_t> chimeric_edges;

        // DFS
        queue.push_front(source);
        visited[visited_length++] = source;
        while (queue.size() != 0 && !found_chim_sink) {
            int32_t v = queue.front();
            queue.pop_front();
            const auto& curr_node = nodes_[v];

            for (const auto& edge: curr_node->suffix_edges) {
                int32_t w = edge->end_node->id;
                if (predecessor[w] != -1 || w == source) {
                    // Cycle or bubble
                    continue;
                }

                visited[visited_length++] = w;
                predecessor[w] = v;
                queue.push_front(w);

                int32_t w_pair = edge->end_node->pair->id;
                if (predecessor[w_pair] != -1) {
                    // Chimeric link!
                    chim_sink = w;
                    chim_sink_pair = w_pair;
                    found_chim_sink = true;
                    break;
                }
            }
        }

        if (found_chim_sink) {
            // fprintf(stderr, "Source = %d, %d, %d\n", source, chim_sink, chim_sink_pair);

            std::vector<int32_t> path;
            extract_path(path, chim_sink);
            if (path.back() != chim_sink_pair) {
                std::vector<int32_t> other_path;
                extract_path(other_path, chim_sink_pair);

                int32_t ancestor_i = -1, ancestor_j = -1;
                for (uint32_t i = 0; i < path.size(); ++i) {
                    for (uint32_t j = 0; j < other_path.size(); ++j) {
                        if (path[i] == other_path[j]) {
                            ancestor_i = i;
                            ancestor_j = j;
                            break;
                        }
                    }
                    if (ancestor_i != -1) {
                        break;
                    }
                }
                if (ancestor_i != -1) path.resize(ancestor_i + 1);
                std::reverse(path.begin(), path.end());
                if (ancestor_j != -1) other_path.resize(ancestor_j + 1);

                for (uint32_t j = 1; j < other_path.size(); ++j) {
                    path.push_back(nodes_[other_path[j]]->pair->id);
                }
            } else {
                std::reverse(path.begin(), path.end());
            }

            // check chimeric patterns
            std::vector<uint32_t> chimeric_edges;
            check_path(chimeric_edges, path);

            if (chimeric_edges.size() > 0) {
                // remove chimeric edges
                // for (const auto& pid: path) fprintf(stderr, "%d -> ", pid);
                // fprintf(stderr, "\n");
                for (const auto& edge_id: chimeric_edges) {
                    // fprintf(stderr, "%d\n", edge_id);
                    // fprintf(stderr, "Removing: %d -> %d\n", edges_[edge_id]->begin_node->id, edges_[edge_id]->end_node->id);
                    edges_[edge_id]->mark = true;
                    edges_[edge_id]->pair->mark = true;
                    marked_edges_.insert(edge_id);
                    marked_edges_.insert(edges_[edge_id]->pair->id);
                }
                ++num_chimeric_edges;
                remove_marked_edges();
            } else {
                // for (const auto& pid: path) fprintf(stderr, "%d -> ", pid);
                // fprintf(stderr, "\n");
            }
        }

        queue.clear();
        for (uint32_t i = 0; i < visited_length; ++i) {
            predecessor[visited[i]] = -1;
        }
        visited_length = 0;
    }

    fprintf(stderr, "    removed %u edges\n", num_chimeric_edges);
    fprintf(stderr, "  }\n");

    return num_chimeric_edges;
}

uint32_t Graph::remove_bubbles() {

    fprintf(stderr, "  Bubble removal {\n");

    std::vector<uint32_t> distance(nodes_.size(), 0);
    std::vector<uint32_t> visited(nodes_.size(), 0);
    uint32_t visited_length = 0;
    std::vector<int32_t> predecessor(nodes_.size(), -1);
    std::deque<uint32_t> queue;

    auto extract_path = [&](std::vector<uint32_t>& dst, int32_t source, int32_t sink) -> void {
        int32_t curr_id = sink;
        while (curr_id != source) {
            dst.push_back(curr_id);
            curr_id = predecessor[curr_id];
        }
        dst.push_back(source);
        std::reverse(dst.begin(), dst.end());
    };

    auto is_valid_bubble = [](const std::vector<uint32_t>& path, const std::vector<uint32_t>& other_path) -> bool {
        std::set<uint32_t> node_set;
        for (const auto& id: path) node_set.insert(id);
        for (const auto& id: other_path) node_set.insert(id);
        if (path.size() + other_path.size() - 2 != node_set.size()) {
            return false;
        }
        for (const auto& id: path) {
            uint32_t pair_id = (id % 2 == 0) ? id + 1 : id - 1;
            if (node_set.count(pair_id) != 0) {
                return false;
            }
        }
        return true;
    };

    auto path_params = [&](const std::vector<uint32_t>& path, uint32_t& num_reads, uint32_t& matching_bases, double& quality) -> void {
        num_reads = 0;
        for (const auto& it: path) num_reads += nodes_[it]->unitig_size;
        for (const auto& edge: nodes_[path[0]]->suffix_edges) {
            if (edge->end_node->id == path[1]) {
                quality = edge->quality;
                matching_bases = edge->matching_bases();
                break;
            }
        }
        return;
    };

    auto remove_path = [&](const std::vector<uint32_t>& path) -> void {
        for (uint32_t i = 0; i < path.size() - 1; ++i) {
            const auto& node = nodes_[path[i]];
            if (i != 0 && (node->in_degree() > 1 || node->out_degree() > 1)) {
                // fprintf(stderr, "Breaking at %d\n", path[i]);
                break;
            }
            for (const auto& edge: node->suffix_edges) {
                if (edge->end_node->id == path[i + 1]) {
                    edge->mark = true;
                    edge->pair->mark = true;
                    marked_edges_.insert(edge->id);
                    marked_edges_.insert(edge->pair->id);
                    break;
                }
            }
        }
        return;
    };

    uint32_t num_bubbles_popped = 0;

    std::vector<uint32_t> bubble_candidates;
    locate_bubble_sources(bubble_candidates);
    for (const auto& id: bubble_candidates) {
        const auto& node = nodes_[id];
        // fprintf(stderr, "Considering bubble source candidate %d for bubble popping\r", id);

        // for (const auto& node: nodes_) {
        if (node == nullptr || node->out_degree() < 2) continue;

        bool found_sink = false;
        int32_t sink = 0, sink_other_predecesor = 0;
        int32_t source = node->id;

        // BFS
        queue.push_back(source);
        visited[visited_length++] = source;
        while (queue.size() != 0 && !found_sink) {
            int32_t v = queue.front();
            const auto& curr_node = nodes_[v];

            queue.pop_front();

            for (const auto& edge: curr_node->suffix_edges) {
                int32_t w = edge->end_node->id;

                if (w == source) {
                    // Cycle
                    continue;
                }

                if (distance[v] + edge->length > kMaxBubbleLength) {
                    // Out of reach
                    continue;
                }

                distance[w] = distance[v] + edge->length;
                visited[visited_length++] = w;
                queue.push_back(w);

                if (predecessor[w] != -1) {
                    sink = w;
                    sink_other_predecesor = v;
                    found_sink = true;
                    break;
                }

                predecessor[w] = v;
            }
        }

        if (found_sink) {
            // fprintf(stderr, "Source = %u, sink = %u, sink_predecesors = [%u, %u]\n", source, sink, predecessor[sink], sink_other_predecesor);

            std::vector<uint32_t> path, other_path;
            extract_path(path, source, sink);
            other_path.push_back(sink);
            extract_path(other_path, source, sink_other_predecesor);

            /*fprintf(stderr, "Path 1:");
            for (const auto& it: path) fprintf(stderr, " %d", it);
            fprintf(stderr, "\n");

            fprintf(stderr, "Path 2:");
            for (const auto& it: other_path) fprintf(stderr, " %d", it);
            fprintf(stderr, "\n");*/

            if (!is_valid_bubble(path, other_path)) {
                // fprintf(stderr, "Not valid bubble!\n");
            } else {
                uint32_t path_reads = 0, path_matching_bases = 0;
                double path_quality = 0;
                path_params(path, path_reads, path_matching_bases, path_quality);

                uint32_t other_path_reads = 0, other_path_matching_bases = 0;
                double other_path_quality = 0;
                path_params(other_path, other_path_reads, other_path_matching_bases,
                    other_path_quality);

                // fprintf(stderr, "Path 1 = (%d, %d, %g) | Path 2 = (%d, %d, %g) | Worse path is ",
                //      path_reads, path_matching_bases, path_quality,
                //      other_path_reads, other_path_matching_bases, other_path_quality);

                if (path_reads > other_path_reads || (path_reads == other_path_reads && path_matching_bases > other_path_matching_bases)) {
                    // fprintf(stderr, "2\n");
                    remove_path(other_path);
                } else {
                    // fprintf(stderr, "1\n");
                    remove_path(path);
                }
                remove_marked_edges();
                ++num_bubbles_popped;
            }
        }

        queue.clear();
        for (uint32_t i = 0; i < visited_length; ++i) {
            distance[visited[i]] = 0;
            predecessor[visited[i]] = -1;
        }
        visited_length = 0;
    }

    remove_isolated_nodes();
    // if (bubble_candidates.size() != 0) fprintf(stderr, "\n");

    fprintf(stderr, "    popped %d bubbles (out of %zu)\n", num_bubbles_popped, bubble_candidates.size());
    fprintf(stderr, "  }\n");

    return num_bubbles_popped;
}

uint32_t Graph::create_unitigs() {

    fprintf(stderr, "  Creating unitigs {\n");

    uint32_t node_id = nodes_.size();
    std::vector<bool> visited(nodes_.size(), false);
    std::vector<std::unique_ptr<Node>> new_nodes;

    uint32_t num_unitigs_created = 0;

    for (const auto& node: nodes_) {
        if (node == nullptr || visited[node->id] || node->is_junction()) continue;
        // fprintf(stderr, "Considering node %d for unittiging\r", node->id);

        bool is_circular = false;
        auto bnode = node.get();
        while (!bnode->is_junction()) {
            visited[bnode->id] = true;
            visited[bnode->pair->id] = true;
            if (bnode->in_degree() == 0 || bnode->prefix_edges.front()->begin_node->is_junction()) {
                break;
            }
            bnode = bnode->prefix_edges.front()->begin_node;
            if (bnode->id == node->id) {
                is_circular = true;
                break;
            }
        }

        auto enode = node.get();
        while (!enode->is_junction()) {
            visited[enode->id] = true;
            visited[enode->pair->id] = true;
            if (enode->out_degree() == 0 || enode->suffix_edges.front()->end_node->is_junction()) {
                break;
            }
            //++unitig_size;
            enode = enode->suffix_edges.front()->end_node;
            if (enode->id == node->id) {
                is_circular = true;
                break;
            }
        }

        // normal
        Node* unitig = nullptr;
        // reverse_complement
        Node* _unitig = nullptr;

        if (is_circular) {
            unitig = new Node(node_id++, bnode, marked_edges_);
            _unitig = new Node(node_id++, bnode->pair, marked_edges_);
        } else if (bnode->id != enode->id) {
            unitig = new Node(node_id++, bnode, enode, marked_edges_);
            _unitig = new Node(node_id++, enode->pair, bnode->pair, marked_edges_);
        }

        if (unitig != nullptr && _unitig != nullptr) {
            unitig->pair = _unitig;
            _unitig->pair = unitig;

            new_nodes.push_back(std::unique_ptr<Node>(unitig));
            new_nodes.push_back(std::unique_ptr<Node>(_unitig));

            // fprintf(stderr, "Unitig: %d -> %d && %d -> %d\n", bnode->id, enode->id, enode->pair->id, bnode->pair->id);
            ++num_unitigs_created;
        }
    }

    for (auto& node: new_nodes) {
        nodes_.push_back(std::move(node));
    }

    remove_marked_edges();
    remove_isolated_nodes();

    fprintf(stderr, "    created %u new unitigs\n", num_unitigs_created);
    fprintf(stderr, "  }\n");

    return num_unitigs_created;
}

void Graph::print_contigs() const {

    fprintf(stderr, "  Contig creation {\n");

    uint32_t contig_id = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr || node->id % 2 == 0 || node->unitig_size < kMinUnitigSize) continue;
        fprintf(stderr, "    >Contig_%d_(Utg:%u), length = %zu (%d -> %d)\n", contig_id,
            node->unitig_size, node->sequence.size(), node->read_ids.front(),
            node->read_ids.back());
        fprintf(stdout, ">Contig_%u_(Utg=%u:Len=%lu)\n%s\n", contig_id++, node->unitig_size, node->sequence.size(), node->sequence.c_str());
    }

    fprintf(stderr, "  }\n");
}

void Graph::locate_bubble_sources(std::vector<uint32_t>& dst) {

    // 0 - unmarked, 1 - in que, 2 - marked
    std::vector<uint8_t> marks(nodes_.size(), false);
    std::deque<uint32_t> node_que;

    for (const auto& node: nodes_) {
        if (node == nullptr || marks[node->id] == true) {
            continue;
        }
        // fprintf(stderr, "Considering node %d as bubble source\r", node->id);

        node_que.push_back(node->id);
        while (!node_que.empty()) {
            const auto& curr_node = nodes_[node_que.front()];
            node_que.pop_front();
            //fprintf(stderr, "%d, %d\n", curr_node->id, marks[curr_node->id]);

            if (marks[curr_node->id] != 2) {
                if (curr_node->out_degree() > 1) {
                    dst.emplace_back(curr_node->id);
                }
                marks[curr_node->id] = 2;
                marks[curr_node->pair->id] = 2;

                for (const auto& edge: curr_node->prefix_edges) {
                    if (marks[edge->begin_node->id] == 0) {
                        marks[edge->begin_node->id] = 1;
                        node_que.push_back(edge->begin_node->id);
                    }
                }
                for (const auto& edge: curr_node->suffix_edges) {
                    if (marks[edge->end_node->id] == 0) {
                        marks[edge->end_node->id] = 1;
                        node_que.push_back(edge->end_node->id);
                    }
                }
            }
        }
    }

    // fprintf(stderr, "\nBubble candidates found (%zu)!\n", dst.size());

    /*for (const auto& it: dst) {
        fprintf(stderr, "%d ", it);
    }
    fprintf(stderr, "\n");*/
}

void Graph::remove_marked_edges() {

    auto delete_edges = [&](std::list<Edge*>& edges) -> void {
        auto edge = edges.begin();
        while (edge != edges.end()) {
            if ((*edge)->mark == true) {
                edge = edges.erase(edge);
            } else {
                ++edge;
            }
        }
    };

    std::unordered_set<uint32_t> marked_nodes;
    for (const auto& it: marked_edges_) {
        marked_nodes.insert(edges_[it]->begin_node->id);
        marked_nodes.insert(edges_[it]->end_node->id);
    }

    for (const auto& it: marked_nodes) {
        delete_edges(nodes_[it]->prefix_edges);
        delete_edges(nodes_[it]->suffix_edges);
    }

    for (const auto& it: marked_edges_) {
        edges_[it].reset();
    }
    marked_edges_.clear();

    /*for (const auto& node: nodes_) {
        if (node == nullptr) continue;
        delete_edges(node->prefix_edges);
        delete_edges(node->suffix_edges);
    }
    for (auto& edge: edges_) {
        if (edge == nullptr) continue;
        if (edge->mark == true) {
            edge.reset();
        }
    }*/
}

void Graph::print_csv(std::string path) const {

    auto graph_file = fopen(path.c_str(), "w");

    for (const auto& node: nodes_) {
        if (node == nullptr || node->id % 2 == 0) continue;
        fprintf(graph_file, "%u [%u] {%d} U:%d,%u [%u] {%d} U:%d,0,-\n",
            node->id, node->length(), node->read_id, node->unitig_size,
            node->pair->id, node->pair->length(), node->pair->read_id, node->pair->unitig_size);
    }

    for (const auto& edge: edges_) {
        if (edge == nullptr) continue;
        fprintf(graph_file, "%u [%u] {%d} U:%d,%u [%u] {%d} U:%d,1,%d %d %g\n",
            edge->begin_node->id, edge->begin_node->length(), edge->begin_node->read_id, edge->begin_node->unitig_size,
            edge->end_node->id, edge->end_node->length(), edge->end_node->read_id, edge->end_node->unitig_size,
            edge->id, edge->length, edge->quality);
    }

    fclose(graph_file);
}

void Graph::remove_selected_nodes_and_edges() {

    uint32_t num_chimeras = 0;
    for (const auto& node: nodes_) {
        if (node == nullptr) continue;
        if (// scerevisiae ont
            node->read_id == 6466 || node->read_id == 4923 || node->read_id == 74242 ||
            node->read_id == 5654 || node->read_id == 58114 || node->read_id == 45102 ||
            node->read_id == 59780 || node->read_id == 9227 || node->read_id == 76003) {

            // ecoli pacbio
            //node->read_id == 41549 || node->read_id == 50611 || node->read_id == 62667 ||
            //node->read_id == 64916 || node->read_id == 65784 || node->read_id == 83480 ||
            //node->read_id == 85123 || node->read_id == 39991 || node->read_id == 1429) {

            node->mark = true;
            node->pair->mark = true;
            for (const auto& edge: node->suffix_edges) {
                edge->mark = true;
                edge->pair->mark = true;
                marked_edges_.insert(edge->id);
                marked_edges_.insert(edge->pair->id);
            }
            for (const auto& edge: node->prefix_edges) {
                edge->mark = true;
                edge->pair->mark = true;
                marked_edges_.insert(edge->id);
                marked_edges_.insert(edge->pair->id);
            }

            ++num_chimeras;
        }
    }

    fprintf(stderr, "Num selected nodes = %u\n", num_chimeras);

    for (const auto& edge: edges_) {
        if (edge == nullptr) continue;
        if (// scerevisiae ont
            edge->id == 79150 || edge->id == 62848 || edge->id == 56114 ||
            edge->id == 44052 || edge->id == 5675 || edge->id == 2616 ||
            edge->id == 33196 || edge->id == 76055 || edge->id == 13810 ||
            edge->id == 67515 || edge->id == 44976 || edge->id == 30040 ||
            edge->id == 50666 || edge->id == 21391 || edge->id == 103954 ||
            edge->id == 118028 || edge->id == 26160 || edge->id == 30568 ||
            edge->id == 57100 || edge->id == 40926 || edge->id == 63575 ||
            edge->id == 47016 || edge->id == 50620 || edge->id == 27678 ||
            edge->id == 38658 || edge->id == 27431) {
            edge->mark = true;
            edge->pair->mark = true;
            marked_edges_.insert(edge->id);
            marked_edges_.insert(edge->pair->id);
        }
    }
    remove_marked_edges();
    remove_isolated_nodes();
}

}
