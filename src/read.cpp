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

constexpr uint32_t kMinReadLength = 1000;
constexpr uint32_t kMinCoverage = 3;
constexpr double kMinDivergence = 0.01;
constexpr double kPitSlopeRatio = 1.817;
constexpr double kHillSlopeRatio = 1.3;
constexpr uint32_t kSlopeWidth = 500;
constexpr double kSlopeWidthRatio = 0.05;
constexpr double kHillWidthRatio = 0.85;
constexpr double kMedianRatio = 2;

using CoverageSubgraph = std::deque<std::pair<int32_t, int32_t>>;

void coverage_subgraph_add(CoverageSubgraph& subgraph, int32_t value,
    int32_t position) {

    while (!subgraph.empty() && subgraph.back().second <= value) {
        subgraph.pop_back();
    }
    subgraph.emplace_back(position, value);
}

void coverage_subgraph_update(CoverageSubgraph& subgraph, int32_t position) {

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
        : id_(id), begin_(0), end_(read_length), coverage_median_(0),
        coverage_graph_(end_ - begin_ + 1, 0), corrected_coverage_graph_(),
        coverage_hills_() {
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

    if (end - begin < kMinReadLength) {
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
        uint32_t divergence = abs(end - begin - other_end + other_begin);

        if (divergence > correction_length * kMinDivergence) {
            continue;
        }

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
        if (!found_begin && coverage_graph_[i] >= kMinCoverage) {
            current_begin = i;
            found_begin = true;
        } else if (found_begin && coverage_graph_[i] < kMinCoverage) {
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

bool ReadInfo::find_coverage_pits(uint16_t dataset_median) {

    CoverageSubgraph left_subgraph, right_subgraph;
    std::vector<int32_t> slopes;

    int32_t k = std::max(kSlopeWidth, static_cast<uint32_t>(kSlopeWidthRatio *
        (end_ - begin_)));
    int32_t read_length = coverage_graph_.size();
    int32_t median_threshold = dataset_median / 2;

    for (int32_t i = -1 * k + 2; i < read_length - 1; ++i) {
        if (i < read_length - k) {
            coverage_subgraph_add(right_subgraph, coverage_graph_[i + k], i + k);
        }
        coverage_subgraph_update(right_subgraph, i);

        if (i > 0) {
            coverage_subgraph_add(left_subgraph, coverage_graph_[i - 1], i - 1);
            coverage_subgraph_update(left_subgraph, i - 1 - k);

            if (coverage_graph_[i] > (uint16_t) median_threshold) {
                continue;
            }

            int32_t current = coverage_graph_[i] * kPitSlopeRatio;
            if (left_subgraph.front().second > current) {
                slopes.push_back(i << 1 | 0);
            }
            if (!right_subgraph.empty() &&
                right_subgraph.front().second > current) {

                slopes.push_back(i << 1 | 1);
            }
        }
    }

    if (slopes.size() > 0) {

        bool is_chimeric = false;
        std::sort(slopes.begin(), slopes.end());

        std::vector<uint32_t> breaking_points(1, begin_);
        for (uint32_t i = 0; i < slopes.size() - 1; ++i) {
            if (!(slopes[i] & 1) && (slopes[i + 1] & 1) &&
                (slopes[i + 1] >> 1) - (slopes[i] >> 1) < k) {

                is_chimeric = true;
                breaking_points.push_back(((slopes[i] >> 1) +
                    (slopes[i + 1] >> 1)) / 2);
            }
        }
        breaking_points.push_back(end_);

        if (is_chimeric) {
            uint32_t new_begin = 0, new_end = 0;
            for (uint32_t i = 0; i < breaking_points.size() - 1; ++i) {
                if (breaking_points[i + 1] - breaking_points[i] >
                    new_end - new_begin) {

                    new_begin = breaking_points[i];
                    new_end = breaking_points[i + 1];
                }
            }
            return shrink_coverage_graph(new_begin, new_end);
        }
    }

    return true;
}

void ReadInfo::find_coverage_hills(uint16_t dataset_median) {

    CoverageSubgraph left_subgraph, right_subgraph;
    std::vector<int32_t> slopes;

    if (!corrected_coverage_graph_.empty()) {
        corrected_coverage_graph_.swap(coverage_graph_);
    }

    int32_t k = std::max(kSlopeWidth, uint32_t(kSlopeWidthRatio *
        (end_ - begin_)));
    int32_t read_length = coverage_graph_.size();

    for (int32_t i = -1 * k + 2; i < read_length; ++i) {
        if (i < read_length - k) {
            coverage_subgraph_add(right_subgraph, coverage_graph_[i + k], i + k);
        }
        coverage_subgraph_update(right_subgraph, i);

        if (i == 0) {
            int32_t current = coverage_graph_[i] * kHillSlopeRatio;
            if (coverage_graph_[i + 1] > dataset_median &&
                !right_subgraph.empty() &&
                right_subgraph.front().second > (int32_t) dataset_median &&
                right_subgraph.front().second > current) {

                slopes.push_back(i << 1 | 1);
            }
        }

        if (i > 0) {
            coverage_subgraph_add(left_subgraph, coverage_graph_[i - 1], i - 1);
            coverage_subgraph_update(left_subgraph, i - 1 - k);

            int32_t current = coverage_graph_[i] * kHillSlopeRatio;
            if (coverage_graph_[i - 1] > dataset_median &&
                left_subgraph.front().second > (int32_t) dataset_median &&
                left_subgraph.front().second > current) {
                slopes.push_back(i << 1 | 0);
            }
            if ((i == read_length - 1 || coverage_graph_[i + 1] > dataset_median) &&
                !right_subgraph.empty() &&
                right_subgraph.front().second > (int32_t) dataset_median &&
                right_subgraph.front().second > current) {
                slopes.push_back(i << 1 | 1);
            }
        }
    }

    if (slopes.size() > 1) {

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
                        slope_regions.emplace_back(fdownslope << 1 | 0,
                            ldownslope);
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

        if (found_fus) {
            slope_regions.emplace_back(fupslope << 1 | 1, lupslope);
        }
        if (found_fds) {
            slope_regions.emplace_back(fdownslope << 1 | 0, ldownslope);
        }
        std::sort(slope_regions.begin(), slope_regions.end());

        // chimeric check - TODO: rerranging slopes does not work all the time
        auto rearrange_slopes = [&](uint32_t i, uint32_t j) -> void {
            uint32_t begin = std::max(slope_regions[i].first >> 1,
                slope_regions[j].first >> 1);
            uint32_t end = std::min(slope_regions[i].second,
                slope_regions[j].second);
            uint32_t min_left_id = begin, min_right_id = begin;
            for (uint32_t s = begin + 1; s < end; ++s) {
                if (coverage_graph_[s] < coverage_graph_[min_left_id]) {
                    min_left_id = s;
                }
                if (coverage_graph_[s] <= coverage_graph_[min_right_id]) {
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
            if ((slope_regions[s].first & 1) &&
                !(slope_regions[s + 1].first & 1) &&
                slope_regions[s].second > (slope_regions[s + 1].first >> 1)) {

                rearrange_slopes(s, s + 1);
            }
        }

        for (uint32_t s = 0; s < slope_regions.size() - 1; ++s) {
            if (slope_regions[s].second > (slope_regions[s + 1].first >> 1)) {
                rearrange_slopes(s, s + 1);
            }
        }

        // find hills TODO: better checker is needed
        auto check_hill = [&](uint32_t begin, uint32_t end, double median) -> bool {
            uint32_t valid_bases = 0;
            for (uint32_t i = begin; i < end; ++i) {
                if (coverage_graph_[i] >= median * kHillSlopeRatio) {
                    ++valid_bases;
                }
            }
            if (valid_bases > 0.85 * (end - begin)) {
                return true;
            }
            return false;
        };

        bool print = false;
        uint32_t max_width = (end_ - begin_) * kHillWidthRatio;
        for (uint32_t r = 0; r < slope_regions.size() - 1; ++r) {
            if ((slope_regions[r].first & 1) &&
                !(slope_regions[r + 1].first & 1)) {
                if (slope_regions[r + 1].second -
                    (slope_regions[r].first >> 1) < max_width) {

                    if (check_hill(slope_regions[r].second + 1,
                        (slope_regions[r + 1].first >> 1) - 1,
                        dataset_median) == false) {

                        continue;
                    }

                    if ((slope_regions[r].first >> 1) < 0.05 * (end_ - begin_) + begin_) {
                        // left hill
                        coverage_hills_.emplace_back(slope_regions[r].first >> 1,
                            slope_regions[r + 1].second);
                    } else if (slope_regions[r + 1].second > 0.95 * (end_ - begin_) + begin_) {
                        // right hill
                        coverage_hills_.emplace_back(slope_regions[r].first >> 1,
                            slope_regions[r + 1].second);
                    }
                }
            }
            if (print) {
                print_csv("graphs/h" + std::to_string(id_), 0);
            }
        }
    }

    if (!corrected_coverage_graph_.empty()) {
        corrected_coverage_graph_.swap(coverage_graph_);
    }
}

bool ReadInfo::is_valid_overlap(uint32_t begin, uint32_t end) const {

    begin += begin_;
    end += begin_;
    uint32_t fuzz = 0.05 * (end_ - begin_);

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

    std::ofstream out(path);
    out << "x " << id_ << " slopes median dataset_median y" << std::endl;
    for (uint32_t i = 0; i < coverage_graph_.size(); ++i) {
        out << i << " " << coverage_graph_[i] << " " <<
            (uint16_t) slope_graph[i] << " " << coverage_median_ << " " <<
            dataset_median << " " << corrected_coverage_graph_[i] << std::endl;
    }
    out.close();
}

}
