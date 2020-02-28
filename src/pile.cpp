/*!
 * @file pile.cpp
 *
 * @brief Pile class source file
 */

#include <algorithm>
#include <sstream>
#include <deque>

#include "overlap.hpp"
#include "pile.hpp"

namespace rala {

using Subpile = std::deque<std::pair<int32_t, int32_t>>;

void subpileAdd(Subpile& src, int32_t value, int32_t position) {
    while (!src.empty() && src.back().second <= value) {
        src.pop_back();
    }
    src.emplace_back(position, value);
}

void subpileUpdate(Subpile& src, int32_t position) {
    while (!src.empty() && src.front().first <= position) {
        src.pop_front();
    }
}

void intervalMerge(std::vector<std::pair<uint32_t, uint32_t>>& intervals) {

    std::vector<std::pair<uint32_t, uint32_t>> tmp;
    std::vector<bool> is_merged(intervals.size(), false);
    for (uint32_t i = 0; i < intervals.size(); ++i) {
        if (is_merged[i]) {
            continue;
        }
        for (uint32_t j = 0; j < intervals.size(); ++j) {
            if (i != j && !is_merged[j] &&
                intervals[i].first < intervals[j].second &&
                intervals[i].second > intervals[j].first) {

                is_merged[j] = true;
                intervals[i].first = std::min(intervals[i].first, intervals[j].first);
                intervals[i].second = std::max(intervals[i].second, intervals[j].second);
            }
        }
        tmp.emplace_back(intervals[i].first, intervals[i].second);
    }
    intervals.swap(tmp);
}

std::unique_ptr<Pile> createPile(uint64_t id, uint32_t read_length) {
    return std::unique_ptr<Pile>(new Pile(id, read_length));
}

Pile::Pile(uint64_t id, uint32_t read_length)
        : id_(id), begin_(0), end_(read_length), p10_(0), median_(0),
        data_(end_ - begin_, 0), repeat_hills_(), repeat_hill_coverage_(),
        chimeric_pits_(), chimeric_hills_(), chimeric_hill_coverage_(), points_() {
}

std::vector<std::pair<uint32_t, uint32_t>> Pile::find_slopes(double q) {

    std::vector<std::pair<uint32_t, uint32_t>> slope_regions;

    int32_t k = 847;
    int32_t read_length = data_.size();

    Subpile left_subpile;
    uint32_t first_down = 0, last_down = 0;
    bool found_down = false;

    Subpile right_subpile;
    uint32_t first_up = 0, last_up = 0;
    bool found_up = false;

    // find slope regions
    for (int32_t i = 0; i < k; ++i) {
        subpileAdd(right_subpile, data_[i], i);
    }
    for (int32_t i = 0; i < read_length; ++i) {
        if (i > 0) {
            subpileAdd(left_subpile, data_[i - 1], i - 1);
        }
        subpileUpdate(left_subpile, i - 1 - k);

        if (i < read_length - k) {
            subpileAdd(right_subpile, data_[i + k], i + k);
        }
        subpileUpdate(right_subpile, i);

        int32_t current_value = data_[i] * q;
        if (i != 0 && left_subpile.front().second > current_value) {
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
        if (i != (read_length - 1) && right_subpile.front().second > current_value) {
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

    if (slope_regions.empty()) {
        return slope_regions;
    }

    while (true) {
        std::sort(slope_regions.begin(), slope_regions.end());

        bool is_changed = false;
        for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
            if (slope_regions[i].second < (slope_regions[i + 1].first >> 1)) {
                continue;
            }

            std::vector<std::pair<uint32_t, uint32_t>> subregions;
            if (slope_regions[i].first & 1) {
                right_subpile.clear();
                found_up = false;

                uint32_t subpile_begin = slope_regions[i].first >> 1;
                uint32_t subpile_end = std::min(slope_regions[i].second,
                    slope_regions[i + 1].second);

                for (uint32_t j = subpile_begin; j < subpile_end + 1; ++j) {
                    subpileAdd(right_subpile, data_[j], j);
                }
                for (uint32_t j = subpile_begin; j < subpile_end; ++j) {
                    subpileUpdate(right_subpile, j);
                    if (data_[j] * q < right_subpile.front().second) {
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
                slope_regions[i].first = subpile_end << 1 | 1;

            } else {
                if (slope_regions[i].second == (slope_regions[i + 1].first >> 1)) {
                    continue;
                }

                left_subpile.clear();
                found_down = false;

                uint32_t subpile_begin = std::max(slope_regions[i].first >> 1,
                    slope_regions[i + 1].first >> 1);
                uint32_t subpile_end = slope_regions[i].second;

                for (uint32_t j = subpile_begin; j < subpile_end + 1; ++j) {
                    if (!left_subpile.empty() && data_[j] * q < left_subpile.front().second) {
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
                    subpileAdd(left_subpile, data_[j], j);
                }
                if (found_down) {
                    subregions.emplace_back(first_down, last_down);
                }

                for (const auto& it: subregions) {
                    slope_regions.emplace_back(it.first << 1 | 0, it.second);
                }
                slope_regions[i].second = subpile_begin;
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

            uint32_t subpile_begin = slope_regions[i].second;
            uint32_t subpile_end = slope_regions[i + 1].first >> 1;

            if (subpile_end - subpile_begin > static_cast<uint32_t>(k)) {
                continue;
            }

            uint16_t max_subpile_coverage = 0;
            for (uint32_t j = subpile_begin + 1; j < subpile_end; ++j) {
                max_subpile_coverage = std::max(max_subpile_coverage, data_[j]);
            }

            uint32_t last_valid_point = slope_regions[i].first >> 1;
            for (uint32_t j = slope_regions[i].first >> 1; j <= subpile_begin; ++j) {
                if (max_subpile_coverage > data_[j] * q) {
                    last_valid_point = j;
                }
            }

            uint32_t first_valid_point = slope_regions[i + 1].second;
            for (uint32_t j = subpile_end; j <= slope_regions[i + 1].second; ++j) {
                if (max_subpile_coverage > data_[j] * q) {
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

void Pile::find_median() {

    std::vector<uint16_t> valid_data(data_.begin() + begin_, data_.begin() + end_);

    std::nth_element(valid_data.begin(), valid_data.begin() + valid_data.size() / 2,
        valid_data.end());
    median_ = valid_data[valid_data.size() / 2];

    std::nth_element(valid_data.begin(), valid_data.begin() + valid_data.size() / 10,
        valid_data.end());
    p10_ = valid_data[valid_data.size() / 10];
}

void Pile::add_layers(std::vector<uint32_t>& overlap_bounds) {

    if (overlap_bounds.empty()) {
        return;
    }

    for (std::uint32_t i = 0; i < overlap_bounds.size(); i += 2) {
        points_.emplace_back(overlap_bounds[i] >> 1, overlap_bounds[i + 1] >> 1);
    }
    std::sort(points_.begin(), points_.end());

    std::sort(overlap_bounds.begin(), overlap_bounds.end());

    uint16_t coverage = 0;
    uint32_t last_bound = begin_;
    for (const auto& bound: overlap_bounds) {
        if (coverage > 0) {
            for (uint32_t i = last_bound; i < (bound >> 1); ++i) {
                data_[i] += coverage;
            }
        }
        last_bound = (bound >> 1);
        if (bound & 1) {
            --coverage;
        } else {
            ++coverage;
        }
    }
}

bool Pile::shrink(uint32_t begin, uint32_t end) {

    if (begin > end) {
        fprintf(stderr, "[rala::Pile::shrink] error: "
            "invalid begin, end coordinates!\n");
        exit(1);
    }

    if (end - begin < 1260) {
        return false;
    }

    for (uint32_t i = begin_; i < begin; ++i) {
        data_[i] = 0;
    }
    begin_ = begin;

    for (uint32_t i = end; i < end_; ++i) {
        data_[i] = 0;
    }
    end_ = end;

    return true;
}

bool Pile::find_valid_region() {

    uint32_t new_begin = 0, new_end = 0, current_begin = 0;
    bool found_begin = false;
    for (uint32_t i = begin_; i < end_; ++i) {
        if (!found_begin && data_[i] >= 4) {
            current_begin = i;
            found_begin = true;
        } else if (found_begin && data_[i] < 4) {
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

    return shrink(new_begin, new_end);
}

void Pile::find_chimeric_pits() {

    auto slope_regions = find_slopes(1.82);
    if (slope_regions.empty()) {
        return;
    }

    for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
        if (!(slope_regions[i].first & 1) && (slope_regions[i + 1].first & 1)) {
            chimeric_pits_.emplace_back(slope_regions[i].first >> 1,
                slope_regions[i + 1].second);
        }
    }
    intervalMerge(chimeric_pits_);
}

bool Pile::break_over_chimeric_pits(uint16_t dataset_median) {

    auto is_chimeric_pit = [&](uint32_t begin, uint32_t end) -> bool {
        for (uint32_t i = begin; i <= end; ++i) {
            if (data_[i] * 1.84 <= dataset_median) {
                return true;
            }
        }
        return false;
    };

    uint32_t begin = 0, end = 0, last_begin = this->begin_;
    std::vector<std::pair<uint32_t, uint32_t>> tmp;

    for (const auto& it: chimeric_pits_) {
        if (begin_ > it.first || end_ < it.second) {
            continue;
        }
        if (is_chimeric_pit(it.first, it.second)) {
            if (it.first - last_begin > end - begin) {
                begin = last_begin;
                end = it.first;
            }
            last_begin = it.second;
        } else {
            tmp.emplace_back(it);
        }
    }
    if (this->end_ - last_begin > end - begin) {
        begin = last_begin;
        end = this->end_;
    }

    chimeric_pits_.swap(tmp);

    return shrink(begin, end);
}

void Pile::find_chimeric_hills() {

    auto slope_regions = find_slopes(1.3);
    if (slope_regions.empty()) {
        return;
    }

    auto is_chimeric_hill = [&](
        const std::pair<uint32_t, uint32_t>& begin,
        const std::pair<uint32_t, uint32_t>& end) -> bool {

        if ((begin.first >> 1) < 0.05 * (this->end_ - this->begin_) + this->begin_ ||
            end.second > 0.95 * (this->end_ - this->begin_) + this->begin_ ||
            (end.first >> 1) - begin.second > 840) {
            return false;
        }

        uint32_t peak_value = 1.3 * std::max(data_[begin.second],
            data_[end.first >> 1]);

        for (uint32_t i = begin.second + 1; i < (end.first >> 1); ++i) {
            if (data_[i] > peak_value) {
                return true;
            }
        }
        return false;
    };

    uint32_t fuzz = 420;
    for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
        if (!(slope_regions[i].first & 1)) {
            continue;
        }

        for (uint32_t j = i + 1; j < slope_regions.size(); ++j) {
            if (slope_regions[j].first & 1) {
                continue;
            }

            if (is_chimeric_hill(slope_regions[i], slope_regions[j])) {
                uint32_t begin = (slope_regions[i].first >> 1) - this->begin_ > fuzz ?
                    (slope_regions[i].first >> 1) - fuzz : this->begin_;
                uint32_t end = this->end_ - slope_regions[j].second > fuzz ?
                    slope_regions[j].second + fuzz : this->end_;
                chimeric_hills_.emplace_back(begin, end);
            }
        }
    }
    intervalMerge(chimeric_hills_);

    chimeric_hill_coverage_.resize(chimeric_hills_.size(), 0);
}

void Pile::check_chimeric_hills(const std::unique_ptr<Overlap>& overlap) {

    uint32_t begin = this->begin_ + (overlap->a_id() == id_ ? overlap->a_begin() :
        overlap->b_begin());
    uint32_t end = this->begin_ + (overlap->a_id() == id_ ? overlap->a_end() :
        overlap->b_end());

    for (uint32_t i = 0; i < chimeric_hills_.size(); ++i) {
        if (begin < chimeric_hills_[i].first && end > chimeric_hills_[i].second) {
            ++chimeric_hill_coverage_[i];
        }
    }
}

bool Pile::break_over_chimeric_hills() {

    uint32_t begin = 0, end = 0, last_begin = this->begin_;

    for (uint32_t i = 0; i < chimeric_hills_.size(); ++i) {
        if (begin_ > chimeric_hills_[i].first || end_ < chimeric_hills_[i].second) {
            continue;
        }
        if (chimeric_hill_coverage_[i] > 3) {
            continue;
        }

        if (chimeric_hills_[i].first - last_begin > end - begin) {
            begin = last_begin;
            end = chimeric_hills_[i].first;
        }
        last_begin = chimeric_hills_[i].second;
    }
    if (this->end_ - last_begin > end - begin) {
        begin = last_begin;
        end = this->end_;
    }

    std::vector<std::pair<uint32_t, uint32_t>>().swap(chimeric_hills_);
    std::vector<uint32_t>().swap(chimeric_hill_coverage_);

    return shrink(begin, end);
}

void Pile::find_repetitive_hills(uint16_t dataset_median) {

    // TODO: remove?
    if (median_ > 1.42 * dataset_median) {
        dataset_median = std::max(dataset_median, p10_);
    }

    auto slope_regions = find_slopes(1.42);
    if (slope_regions.empty()) {
        return;
    }

    auto is_repeat_hill = [&](
        const std::pair<uint32_t, uint32_t>& begin,
        const std::pair<uint32_t, uint32_t>& end) -> bool {

        if (((end.first >> 1) + end.second) / 2 -
            ((begin.first >> 1) + begin.second) / 2 > 0.84 * (this->end_ - this->begin_)) {
            return false;
        }
        bool found_peak = false;
        uint32_t peak_value = 1.42 * std::max(data_[begin.second], data_[end.first >> 1]);
        uint32_t valid_points = 0;
        uint32_t min_value = dataset_median * 1.42;

        for (uint32_t i = begin.second + 1; i < (end.first >> 1); ++i) {
            if (data_[i] > min_value) {
                ++valid_points;
            }
            if (data_[i] > peak_value) {
                found_peak = true;
            }
        }

        if (!found_peak || valid_points < 0.9 * ((end.first >> 1) - begin.second)) {
            return false;
        }
        return true;
    };

    for (uint32_t i = 0; i < slope_regions.size() - 1; ++i) {
        if (!(slope_regions[i].first & 1)) {
            continue;
        }
        for (uint32_t j = i + 1; j < slope_regions.size(); ++j) {
            if (slope_regions[j].first & 1) {
                continue;
            }

            if (is_repeat_hill(slope_regions[i], slope_regions[j])) {
                repeat_hills_.emplace_back(
                    slope_regions[i].second - 0.336 *
                        (slope_regions[i].second - (slope_regions[i].first >> 1)),
                    (slope_regions[j].first >> 1) + 0.336 *
                        (slope_regions[j].second - (slope_regions[j].first >> 1)));
            }
        }
    }

    intervalMerge(repeat_hills_);
    for (auto& it: repeat_hills_) {
        it.first = std::max(begin_, it.first);
        it.second = std::min(end_, it.second);
    }

    repeat_hill_coverage_.resize(repeat_hills_.size(), false);
}

void Pile::check_repetitive_hills(const std::unique_ptr<Overlap>& overlap) {

    uint32_t begin = overlap->b_begin();
    uint32_t end = overlap->b_end();
    uint32_t fuzz = 420;

    for (uint32_t i = 0; i < repeat_hills_.size(); ++i) {
        if (begin < repeat_hills_[i].second && repeat_hills_[i].first < end) {
            if (repeat_hills_[i].first < 0.1 * (this->end_ - this->begin_) + this->begin_ &&
                begin - this->begin_ < this->end_ - end) {
                // left hill
                if (end >= repeat_hills_[i].second + fuzz) {
                    repeat_hill_coverage_[i] = true;
                }
            } else if (repeat_hills_[i].second > 0.9 * (this->end_ - this->begin_) + this->begin_ &&
                begin - this->begin_ > this->end_ - end) {
                // right hill
                if (begin + fuzz <= repeat_hills_[i].first) {
                    repeat_hill_coverage_[i] = true;
                }
            }

        }
    }
}

void Pile::add_repetitive_region(uint32_t begin, uint32_t end) {

    if (begin > data_.size() || end > data_.size()) {
        fprintf(stderr, "[rala::Pile::add_repetitive_region] error: "
            "[begin,end] out of bounds!\n");
        exit(1);
    }

    repeat_hills_.emplace_back(begin, end);
}

bool Pile::is_valid_overlap(uint32_t begin, uint32_t end) const {

    uint32_t fuzz = 420;

    auto check_hills = [&](const std::vector<std::pair<uint32_t, uint32_t>>& hills) -> bool {
        for (uint32_t i = 0; i < hills.size(); ++i) {
            const auto& it = hills[i];
            if (begin < it.second && it.first < end) {
                if (it.first < 0.1 * (this->end_ - this->begin_) + this->begin_) {
                    // left hill
                    if (end < it.second + fuzz && repeat_hill_coverage_[i]) {
                        return false;
                    }
                } else if (it.second > 0.9 * (this->end_ - this->begin_) + this->begin_) {
                    // right hill
                    if (begin + fuzz > it.first && repeat_hill_coverage_[i]) {
                        return false;
                    }
                }
            }
        }
        return true;
    };

    return check_hills(repeat_hills_);
}

std::string Pile::to_json() const {

    std::stringstream ss;
    ss << "\"" << id_ << "\":{";

    ss << "\"y\":[";
    //for (uint32_t i = 0; i < data_.size(); ++i) {
    for (uint32_t i = begin_; i < end_; ++i) {
        ss << data_[i];
        if (i < end_ - 1) {
            ss << ",";
        }
    }
    ss << "],";

    ss << "\"b\":" << begin_ << ",";
    ss << "\"e\":" << end_ << ",";

    /*
    ss << "\"h\":[";
    for (uint32_t i = 0; i < repeat_hills_.size(); ++i) {
        ss << repeat_hills_[i].first << "," << repeat_hills_[i].second;
        if (i < repeat_hills_.size() - 1) {
            ss << ",";
        }
    }
    ss << "],";
    */

    ss << "\"m\":" << median_ << ",";
    ss << "\"p10\":" << p10_ << ",";

    ss << "\"p\":[";
    for (uint32_t i = 0; i < points_.size(); ++i) {
        ss << points_[i].first << "," << points_[i].second;
        if (i < points_.size() - 1) {
            ss << ",";
        }
    }
    ss << "]";
    ss << "}";

    return ss.str();
}

}
