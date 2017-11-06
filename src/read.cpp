/*!
 * @file read.cpp
 *
 * @brief Read class source file
 */

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <deque>

#include "read.hpp"

namespace rala {

enum class CoverageHillType {
    kInvalid,
    kNormal,
    kRepeat,
    kChimeric
};

CoverageHillType checkCoverageHill(const std::pair<uint32_t, uint32_t>& hill_begin,
    const std::pair<uint32_t, uint32_t>& hill_end, uint32_t dataset_median,
    const std::vector<uint16_t>& coverage_graph, uint32_t begin, uint32_t end) {

    if (hill_end.second - (hill_begin.first >> 1) > 0.84 * (end - begin)) {
        return CoverageHillType::kInvalid;
    }

    bool found_peak = false;
    uint32_t peak_value = 1.3 * std::max(coverage_graph[hill_begin.second],
        coverage_graph[hill_end.first >> 1]);

    uint32_t valid_points = 0;
    uint32_t min_value = dataset_median * 1.3;

    for (uint32_t i = hill_begin.second + 1; i < (hill_end.first >> 1); ++i) {
        if (coverage_graph[i] > min_value) {
            ++valid_points;
        }
        if (coverage_graph[i] > peak_value) {
            found_peak = true;
        }
    }
    if (!found_peak) {
        return CoverageHillType::kInvalid;
    }
    if (valid_points < 0.84 * ((hill_end.first >> 1) - hill_begin.second)) {
        if ((hill_end.first >> 1) - hill_begin.second < 596 &&
            (hill_begin.first >> 1) > 0.05 * (end - begin) + begin &&
            (hill_end.second) < 0.95 * (end - begin) + begin) {

            return CoverageHillType::kChimeric;
        }
        return CoverageHillType::kNormal;
    }
    return CoverageHillType::kRepeat;
}

void mergeCoverageHills(std::vector<std::pair<uint32_t, uint32_t>>& hills) {

    std::vector<std::pair<uint32_t, uint32_t>> tmp;
    std::vector<bool> is_merged(hills.size(), false);
    for (uint32_t i = 0; i < hills.size(); ++i) {
        if (is_merged[i]) {
            continue;
        }
        for (uint32_t j = 0; j < hills.size(); ++j) {
            if (i != j && !is_merged[j] &&
                hills[i].first < hills[j].second &&
                hills[i].second > hills[j].first) {

                is_merged[j] = true;
                hills[i].first = std::min(hills[i].first, hills[j].first);
                hills[i].second = std::max(hills[i].second, hills[j].second);
            }
        }
        tmp.emplace_back(hills[i].first, hills[i].second);
    }
    hills.swap(tmp);
};

using CoverageSubgraph = std::deque<std::pair<int32_t, int32_t>>;

void coverageSubgraphAdd(CoverageSubgraph& subgraph, int32_t value,
    int32_t position) {

    while (!subgraph.empty() && subgraph.back().second <= value) {
        subgraph.pop_back();
    }
    subgraph.emplace_back(position, value);
}

void coverageSubgraphUpdate(CoverageSubgraph& subgraph, int32_t position) {
    while (!subgraph.empty() && subgraph.front().first <= position) {
        subgraph.pop_front();
    }
}

Read::Read(uint64_t id, const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length)
        : id_(id), name_(name, name_length), sequence_(sequence,
        sequence_length), reverse_complement_() {
}

Read::Read(uint64_t id, const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length)
        : Read(id, name, name_length, sequence, sequence_length) {
}

void Read::update(uint32_t begin, uint32_t end) {
    sequence_ = sequence_.substr(begin, end - begin);
    if (!reverse_complement_.empty()) {
        create_reverse_complement();
    }
}

void Read::create_reverse_complement() {
    reverse_complement_.clear();
    for (int32_t i = sequence_.size() - 1; i >= 0; --i) {
        switch (sequence_[i]) {
            case 'A':
                reverse_complement_ += 'T';
                break;
            case 'T':
                reverse_complement_ += 'A';
                break;
            case 'C':
                reverse_complement_ += 'G';
                break;
            case 'G':
                reverse_complement_ += 'C';
                break;
            default:
                break;
        }
    }

    assert(reverse_complement_.size() == sequence_.size());
}

std::unique_ptr<ReadInfo> createReadInfo(uint64_t id, uint32_t read_length) {
    return std::unique_ptr<ReadInfo>(new ReadInfo(id, read_length));
}

ReadInfo::ReadInfo(uint64_t id, uint32_t read_length)
        : id_(id), begin_(0), end_(read_length), coverage_p10_(0), coverage_median_(0),
        coverage_graph_(end_ - begin_ + 1, 0), corrected_coverage_graph_(),
        coverage_hills_() {
}

std::vector<std::pair<uint32_t, uint32_t>> ReadInfo::find_coverage_slopes(double q) {

    std::vector<std::pair<uint32_t, uint32_t>> slope_regions;

    int32_t k = 847;
    int32_t read_length = coverage_graph_.size();

    CoverageSubgraph left_subgraph;
    uint32_t first_down = 0, last_down = 0;
    bool found_down = false;

    CoverageSubgraph right_subgraph;
    uint32_t first_up = 0, last_up = 0;
    bool found_up = false;

    // find slope regions
    for (int32_t i = 0; i < k; ++i) {
        coverageSubgraphAdd(right_subgraph, coverage_graph_[i], i);
    }
    for (int32_t i = 0; i < read_length; ++i) {
        if (i > 0) {
            coverageSubgraphAdd(left_subgraph, coverage_graph_[i - 1], i - 1);
        }
        coverageSubgraphUpdate(left_subgraph, i - 1 - k);

        if (i < read_length - k) {
            coverageSubgraphAdd(right_subgraph, coverage_graph_[i + k], i + k);
        }
        coverageSubgraphUpdate(right_subgraph, i);

        int32_t current_value = coverage_graph_[i] * q;
        if (i != 0 && left_subgraph.front().second > current_value) {
            if (found_down) {
                if (i - last_down > 1) {
                    slope_regions.emplace_back(first_down << 1 | 0, last_down);
                    first_down = i;
                }
            } else {
                found_down = true;
                first_down = i;
            }
            last_down = i;
        }
        if (i != (read_length - 1) && right_subgraph.front().second > current_value) {
            if (found_up) {
                if (i - last_up > 1) {
                    slope_regions.emplace_back(first_up << 1 | 1, last_up);
                    first_up = i;
                }
            } else {
                found_up = true;
                first_up = i;
            }
            last_up = i;
        }
    }
    if (found_down) {
        slope_regions.emplace_back(first_down << 1 | 0, last_down);
    }
    if (found_up) {
        slope_regions.emplace_back(first_up << 1 | 1, last_up);
    }

    // rearrange overlapping regions
    while (true) {
        std::sort(slope_regions.begin(), slope_regions.end());

        uint32_t is_changed = false;
        for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
            if (slope_regions[i].second <= (slope_regions[i + 1].first >> 1)) {
                continue;
            }

            std::vector<std::pair<uint32_t, uint32_t>> subregions;
            if (slope_regions[i].first & 1) {
                right_subgraph.clear();
                found_up = false;

                uint32_t subgraph_begin = slope_regions[i].first >> 1;
                uint32_t subgraph_end = std::min(slope_regions[i].second,
                    slope_regions[i + 1].second);

                for (uint32_t j = subgraph_begin; j < subgraph_end + 1; ++j) {
                    coverageSubgraphAdd(right_subgraph, coverage_graph_[j], j);
                }
                for (uint32_t j = subgraph_begin; j < subgraph_end; ++j) {
                    coverageSubgraphUpdate(right_subgraph, j);
                    if (coverage_graph_[j] * q < right_subgraph.front().second) {
                        if (found_up) {
                            if (j - last_up > 1) {
                                subregions.emplace_back(first_up, last_up);
                                first_up = j;
                            }
                        } else {
                            found_up = true;
                            first_up = j;
                        }
                        last_up = j;
                    }
                }
                if (found_up) {
                    subregions.emplace_back(first_up, last_up);
                }

                for (const auto& it: subregions) {
                    slope_regions.emplace_back(it.first << 1 | 1, it.second);
                }
                slope_regions[i].first = subgraph_end << 1 | 1;

            } else {
                left_subgraph.clear();
                found_down = false;

                uint32_t subgraph_begin = std::max(slope_regions[i].first >> 1,
                    slope_regions[i + 1].first >> 1);
                uint32_t subgraph_end = slope_regions[i].second;

                for (uint32_t j = subgraph_begin; j < subgraph_end + 1; ++j) {
                    if (!left_subgraph.empty() &&
                        coverage_graph_[j] * q < left_subgraph.front().second) {

                        if (found_down) {
                            if (j - last_down > 1) {
                                subregions.emplace_back(first_down, last_down);
                                first_down = j;
                            }
                        } else {
                            found_down = true;
                            first_down = j;
                        }
                        last_down = j;
                    }
                    coverageSubgraphAdd(left_subgraph, coverage_graph_[j], j);
                }
                if (found_down) {
                    subregions.emplace_back(first_down, last_down);
                }

                for (const auto& it: subregions) {
                    slope_regions.emplace_back(it.first << 1 | 0, it.second);
                }
                slope_regions[i].second = subgraph_begin;
            }

            is_changed = true;
            break;
        }

        if (!is_changed) {
            break;
        }
    }

    // narrow slope regions
    for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
        if ((slope_regions[i].first & 1) && !(slope_regions[i + 1].first & 1)) {

            uint32_t subgraph_begin = slope_regions[i].second;
            uint32_t subgraph_end = slope_regions[i + 1].first >> 1;

            if (subgraph_end - subgraph_begin > static_cast<uint32_t>(k)) {
                continue;
            }

            uint16_t max_subgraph_coverage = 0;
            for (uint32_t j = subgraph_begin + 1; j < subgraph_end; ++j) {
                max_subgraph_coverage = std::max(max_subgraph_coverage,
                    coverage_graph_[j]);
            }

            uint32_t last_valid_point = slope_regions[i].first >> 1;
            for (uint32_t j = slope_regions[i].first >> 1; j <= subgraph_begin; ++j) {
                if (max_subgraph_coverage > coverage_graph_[j] * q) {
                    last_valid_point = j;
                }
            }

            uint32_t first_valid_point = slope_regions[i + 1].second;
            for (uint32_t j = subgraph_end; j <= slope_regions[i + 1].second; ++j) {
                if (max_subgraph_coverage > coverage_graph_[j] * q) {
                    first_valid_point = j;
                    break;
                }
            }

            slope_regions[i].second = last_valid_point;
            slope_regions[i + 1].first = first_valid_point << 1 | 0;
        }
    }

    return slope_regions;
}

void ReadInfo::update_coverage_graph(std::vector<uint32_t>& shrunken_overlaps) {

    std::sort(shrunken_overlaps.begin(), shrunken_overlaps.end());

    uint16_t coverage = 0;
    uint32_t last_position = 0;
    for (const auto& position: shrunken_overlaps) {
        if (coverage > 0) {
            for (uint32_t i = last_position; i < (position >> 1); ++i) {
                coverage_graph_[i] += coverage;
            }
        }
        last_position = position >> 1;
        if (position & 1) {
            --coverage;
        } else {
            ++coverage;
        }
    }
}

bool ReadInfo::shrink_coverage_graph(uint32_t begin, uint32_t end) {

    assert(begin <= end);

    if (end - begin < 1000) {
        return false;
    }

    for (uint32_t i = begin_; i < begin; ++i) {
        coverage_graph_[i] = 0;
    }
    begin_ = begin;

    for (uint32_t i = end; i < end_; ++i) {
        coverage_graph_[i] = 0;
    }
    end_ = end;

    return true;
}

void ReadInfo::correct_coverage_graph(
    const std::vector<uint32_t>& shrunken_overlaps,
    const std::vector<std::unique_ptr<ReadInfo>>& read_infos) {

    if (corrected_coverage_graph_.empty()) {
        corrected_coverage_graph_ = coverage_graph_;
    }

    for (uint32_t i = 0; i < shrunken_overlaps.size(); i += 6) {
        const auto& other = read_infos[shrunken_overlaps[i + 2]];

        uint32_t begin = begin_ + shrunken_overlaps[i];
        uint32_t end = begin_ + shrunken_overlaps[i + 1];
        uint32_t other_begin = other->begin_ + shrunken_overlaps[i + 3];
        uint32_t other_end = other->begin_ + shrunken_overlaps[i + 4];
        uint32_t orientation = shrunken_overlaps[i + 5];

        assert(other != nullptr);
        assert(begin_ <= begin && begin < end_);
        assert(begin_ < end && end <= end_);
        assert(other->begin_ <= other_begin && other_begin < other->end_);
        assert(other->begin_ < other_end && other_end <= other->end_);

        uint32_t correction_length = std::min(end - begin, other_end -
            other_begin);

        for (uint32_t j = 0; j < correction_length; ++j) {
            if (orientation == 0) {
                corrected_coverage_graph_[begin + j] = std::max(
                    corrected_coverage_graph_[begin + j],
                    other->coverage_graph_[other_begin + j]);
            } else {
                corrected_coverage_graph_[begin + j] = std::max(
                    corrected_coverage_graph_[begin + j],
                    other->coverage_graph_[other_end - j - 1]);
            }
        }
    }
}

void ReadInfo::find_coverage_median() {

    if (!corrected_coverage_graph_.empty()) {
        corrected_coverage_graph_.swap(coverage_graph_);
    }

    std::vector<uint16_t> tmp(coverage_graph_.begin() + begin_,
        coverage_graph_.begin() + end_);
    std::sort(tmp.begin(), tmp.end());

    coverage_p10_ = tmp[tmp.size() / 10];

    coverage_median_ = tmp.size() % 2 == 1 ? tmp[tmp.size() / 2] :
        (tmp[tmp.size() / 2 - 1] + tmp[tmp.size() / 2]) / 2;

    if (!corrected_coverage_graph_.empty()) {
        corrected_coverage_graph_.swap(coverage_graph_);
    }
}

bool ReadInfo::find_valid_region() {

    uint32_t new_begin = 0, new_end = 0, current_begin = 0;
    bool found_begin = false;
    for (uint32_t i = begin_; i < end_; ++i) {
        if (!found_begin && coverage_graph_[i] >= 3) {
            current_begin = i;
            found_begin = true;
        } else if (found_begin && coverage_graph_[i] < 3) {
            if (i - current_begin > new_end - new_begin) {
                new_begin = current_begin;
                new_end = i;
            }
            found_begin = false;
        }
    }
    if (found_begin) {
        if (end_ - current_begin > new_end - new_begin) {
            new_begin = current_begin;
            new_end = end_;
        }
    }

    return this->shrink_coverage_graph(new_begin, new_end);
}

bool ReadInfo::find_chimeric_region(uint16_t dataset_median) {

    if (coverage_median_ > 1.42 * dataset_median) {
        dataset_median = std::max(dataset_median, coverage_p10_);
    }

    // look for chimeric pits
    auto slope_regions = this->find_coverage_slopes(1.817);

    if (!slope_regions.empty()) {
        auto is_chimeric_slope_region = [&](uint32_t begin, uint32_t end) -> bool {
            for (uint32_t i = begin; i < end; ++i) {
                if (coverage_graph_[i] <= dataset_median / 2) {
                    return true;
                }
            }
            return false;
        };

        uint32_t new_begin = 0, new_end = 0, last_slope = begin_;
        for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
            if (!(slope_regions[i].first & 1) && (slope_regions[i + 1].first & 1)) {
                bool is_chimeric =
                    (is_chimeric_slope_region(
                        slope_regions[i].first >> 1,
                        slope_regions[i].second + 1)) |
                    (is_chimeric_slope_region(
                        slope_regions[i + 1].first >> 1,
                        slope_regions[i + 1].second + 1));

                if (is_chimeric) {
                    if ((slope_regions[i].first >> 1) - last_slope >
                        new_end - new_begin) {

                        new_begin = last_slope;
                        new_end = slope_regions[i].first >> 1;
                    }
                    last_slope = slope_regions[i + 1].second;
                }
            }
        }
        if (end_ - last_slope > new_end - new_begin) {
            new_begin = last_slope;
            new_end = end_;
        }
        if (!this->shrink_coverage_graph(new_begin, new_end)) {
            return false;
        }
    }

    // look for chimeric hills
    slope_regions = this->find_coverage_slopes(1.3);

    std::vector<std::pair<uint32_t, uint32_t>> chimeric_hills;
    for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
        if (!(slope_regions[i].first & 1)) {
            continue;
        }
        bool found_hill = false;
        for (uint32_t j = i + 1; j < slope_regions.size(); ++j) {
            if (slope_regions[j].first & 1) {
                if (found_hill) {
                    break;
                }
                continue;
            }
            found_hill = true;
            auto hill_type = checkCoverageHill(slope_regions[i], slope_regions[j],
                dataset_median, coverage_graph_, begin_, end_);

            if (hill_type == CoverageHillType::kChimeric) {
                uint32_t begin = slope_regions[i].second - 0.126 *
                    (slope_regions[i].second - (slope_regions[i].first >> 1));
                uint32_t end = (slope_regions[j].first >> 1) + 0.126 *
                    (slope_regions[j].second - (slope_regions[j].first >> 1));

                chimeric_hills.emplace_back(std::max(begin_, begin),
                    std::min(end_, end));
            }
        }
    }

    if (!chimeric_hills.empty()) {
        mergeCoverageHills(chimeric_hills);

        uint32_t new_begin = begin_, new_end = chimeric_hills.front().first;
        for (uint32_t i = 0; i < chimeric_hills.size() - 1; ++i) {
            if (chimeric_hills[i + 1].first - chimeric_hills[i].second >
                new_end - new_begin) {

                new_begin = chimeric_hills[i].second;
                new_end = chimeric_hills[i + 1].first;
            }
        }
        if (end_ - chimeric_hills.back().second > new_end - new_begin) {
            new_begin = chimeric_hills.back().second;
            new_end = end_;
        }

        return this->shrink_coverage_graph(new_begin, new_end);
    }

    return true;
}

void ReadInfo::find_repetitive_region(uint16_t dataset_median) {

    if (!corrected_coverage_graph_.empty()) {
        corrected_coverage_graph_.swap(coverage_graph_);
    }

    if (coverage_median_ > 1.42 * dataset_median) {
        dataset_median = std::max(dataset_median, coverage_p10_);
    }

    auto slope_regions = this->find_coverage_slopes(1.3);

    std::vector<std::pair<uint32_t, uint32_t>> repeat_hills;
    for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
        if (!(slope_regions[i].first & 1)) {
            continue;
        }
        bool found_hill = false;
        for (uint32_t j = i + 1; j < slope_regions.size(); ++j) {
            if (slope_regions[j].first & 1) {
                if (found_hill) {
                    break;
                }
                continue;
            }
            found_hill = true;
            auto hill_type = checkCoverageHill(slope_regions[i], slope_regions[j],
                dataset_median, coverage_graph_, begin_, end_);

            if (hill_type == CoverageHillType::kRepeat) {
                uint32_t begin = slope_regions[i].second - 0.336 *
                    (slope_regions[i].second - (slope_regions[i].first >> 1));
                uint32_t end = (slope_regions[j].first >> 1) + 0.336 *
                    (slope_regions[j].second - (slope_regions[j].first >> 1));

                repeat_hills.emplace_back(begin, end);
            }
        }
    }

    // store repetitive hills
    mergeCoverageHills(repeat_hills);
    for (const auto& it: repeat_hills) {
        uint32_t hill_begin = std::max(begin_, it.first),
            hill_end = std::min(end_, it.second);
        if (hill_begin < hill_end) {
            coverage_hills_.emplace_back(hill_begin, hill_end);
        }
    }

    if (!corrected_coverage_graph_.empty()) {
        corrected_coverage_graph_.swap(coverage_graph_);
    }
}

bool ReadInfo::is_valid_overlap(uint32_t begin, uint32_t end) const {

    begin += begin_;
    end += begin_;
    uint32_t fuzz = 0.042 * (end_ - begin_);

    for (const auto& it: coverage_hills_) {
        if (begin < it.second && it.first < end) {
            if (it.first < 0.1 * (end_ - begin_) + begin_) {
                // left hill
                if (end < it.second + fuzz) {
                    return false;
                }
            } else if (it.second > 0.9 * (end_ - begin_) + begin_) {
                // right hill
                if (begin > it.first - fuzz) {
                    return false;
                }
            }
        }
    }

    return true;
}

void ReadInfo::print_csv(std::string path, uint16_t dataset_median) const {

    std::vector<uint8_t> slope_graph(coverage_graph_.size(), 0);
    for (const auto& it: coverage_hills_) {
        slope_graph[it.first] = 2;
        slope_graph[it.second] = 1;
    }

    const std::vector<uint16_t>& corrected_coverage_graph =
        corrected_coverage_graph_.empty() ? coverage_graph_ :
        corrected_coverage_graph_;

    std::ofstream out(path);
    out << "x " << id_ << " slopes median dataset_median y" << std::endl;
    for (uint32_t i = 0; i < coverage_graph_.size(); ++i) {
        out << i << " " << coverage_graph_[i] << " " <<
            (uint16_t) slope_graph[i] << " " << coverage_median_ << " " <<
            dataset_median << " " << corrected_coverage_graph[i] << std::endl;
    }
    out.close();
}

}
