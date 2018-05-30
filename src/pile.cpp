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

using Hill = std::pair<uint32_t, uint32_t>;

enum class HillType {
    kInvalid,
    kNormal,
    kRepeat,
    kChimeric
};

HillType hillType(const std::pair<uint32_t, uint32_t>& hill_begin,
    const std::pair<uint32_t, uint32_t>& hill_end, uint32_t dataset_median,
    const std::vector<uint16_t>& pile, uint32_t begin, uint32_t end) {

    if (((hill_end.first >> 1) + hill_end.second) / 2 -
        ((hill_begin.first >> 1) + hill_begin.second) / 2 > 0.84 * (end - begin)) {
        return HillType::kInvalid;
    }

    bool found_peak = false;
    uint32_t peak_value = 1.3 * std::max(pile[hill_begin.second], pile[hill_end.first >> 1]);

    uint32_t valid_points = 0;
    uint32_t min_value = dataset_median * 1.3;

    for (uint32_t i = hill_begin.second + 1; i < (hill_end.first >> 1); ++i) {
        if (pile[i] > min_value) {
            ++valid_points;
        }
        if (pile[i] > peak_value) {
            found_peak = true;
        }
    }
    if (!found_peak) {
        return HillType::kInvalid;
    }
    if (valid_points < 0.84 * ((hill_end.first >> 1) - hill_begin.second)) {
        if ((hill_end.first >> 1) - hill_begin.second < 756 &&
            (hill_begin.first >> 1) > 0.05 * (end - begin) + begin &&
            (hill_end.second) < 0.95 * (end - begin) + begin) {

            return HillType::kChimeric;
        }
        return HillType::kNormal;
    }
    return HillType::kRepeat;
}

void hillMerge(std::vector<std::pair<uint32_t, uint32_t>>& hills) {

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
}

std::unique_ptr<Pile> createPile(uint64_t id, uint32_t read_length) {
    return std::unique_ptr<Pile>(new Pile(id, read_length));
}

Pile::Pile(uint64_t id, uint32_t read_length)
        : id_(id), begin_(0), end_(read_length), p10_(0), median_(0),
        data_(end_ - begin_, 0), corrected_data_(), hills_() {
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

    if (!corrected_data_.empty()) {
        corrected_data_.swap(data_);
    }

    std::vector<uint16_t> valid_data(data_.begin() + begin_, data_.begin() + end_);

    std::nth_element(valid_data.begin(), valid_data.begin() + valid_data.size() / 2,
        valid_data.end());
    median_ = valid_data[valid_data.size() / 2];

    std::nth_element(valid_data.begin(), valid_data.begin() + valid_data.size() / 10,
        valid_data.end());
    p10_ = valid_data[valid_data.size() / 10];

    if (!corrected_data_.empty()) {
        corrected_data_.swap(data_);
    }
}

void Pile::add_layers(std::vector<uint32_t>& overlap_bounds) {

    if (overlap_bounds.empty()) {
        return;
    }

    std::sort(overlap_bounds.begin(), overlap_bounds.end());

    uint16_t coverage = 0;
    uint32_t last_bound = 0;
    for (const auto& bound: overlap_bounds) {
        if (coverage > 0) {
            for (uint32_t i = last_bound; i < (bound >> 1); ++i) {
                data_[i] += coverage;
            }
        }
        last_bound = bound >> 1;
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

    if (end - begin < 1000) {
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

void Pile::correct(const std::vector<std::shared_ptr<Overlap>>& overlaps,
    const std::vector<std::unique_ptr<Pile>>& piles) {

    if (overlaps.empty()) {
        return;
    }

    if (corrected_data_.empty()) {
        corrected_data_ = data_;
    }

    for (const auto& it: overlaps) {

        const auto& other = piles[(it->a_id() == id_ ? it->b_id() : it->a_id())];

        if (other == nullptr) {
            fprintf(stderr, "[rala::Pile::correct] error: missing other pile!\n");
            exit(1);
        }

        uint32_t begin, end, other_begin, other_end;

        if (it->a_id() == id_) {
            begin = begin_ + it->a_begin();
            end = begin_ + it->a_end();
            other_begin = other->begin_ + it->b_begin();
            other_end = other->begin_ + it->b_end();
        } else {
            begin = begin_ + it->b_begin();
            end = begin_ + it->b_end();
            other_begin = other->begin_ + it->a_begin();
            other_end = other->begin_ + it->a_end();
        }

        if (begin < begin_ || begin >= end_ || end <= begin_ || end > end_) {
            fprintf(stderr, "[rala::Pile::correct] error: "
                "invalid begin, end coordinates!\n");
            exit(1);
        }
        if (other_begin < other->begin_ || other_begin >= other->end_ ||
            other_end <= other->begin_ || other_end > other->end_) {
            fprintf(stderr, "[rala::Pile::correct] error: "
                "invalid other begin, end coordinates!\n");
            exit(1);
        }

        uint32_t correction_length = std::min(other_end - other_begin,
            end - begin);

        for (uint32_t i = 0; i < correction_length; ++i) {
            if (it->orientation() == 0) {
                corrected_data_[begin + i] = std::max(corrected_data_[begin + i],
                    other->data_[other_begin + i]);
            } else {
                corrected_data_[begin + i] = std::max(corrected_data_[begin + i],
                    other->data_[other_end - i - 1]);
            }
        }
    }
}

bool Pile::find_valid_region() {

    uint32_t new_begin = 0, new_end = 0, current_begin = 0;
    bool found_begin = false;
    for (uint32_t i = begin_; i < end_; ++i) {
        if (!found_begin && data_[i] >= 3) {
            current_begin = i;
            found_begin = true;
        } else if (found_begin && data_[i] < 3) {
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

bool Pile::find_chimeric_regions(uint16_t dataset_median) {

    // TODO: iterative chimeric detection?
    if (median_ > 1.42 * dataset_median) {
        dataset_median = std::max(dataset_median, p10_);
    }

    // look for chimeric pits
    auto slope_regions = find_slopes(1.817);

    if (!slope_regions.empty()) {
        auto is_chimeric_slope_region = [&](uint32_t begin, uint32_t end) -> bool {
            for (uint32_t i = begin; i < end; ++i) {
                if (data_[i] <= dataset_median / 2) {
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
        if (!shrink(new_begin, new_end)) {
            return false;
        }
    }

    // look for chimeric hills
    slope_regions = find_slopes(1.3);
    if (slope_regions.empty()) {
        return true;
    }

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
            auto hill_type = hillType(slope_regions[i], slope_regions[j],
                dataset_median, data_, begin_, end_);

            if (hill_type == HillType::kChimeric) {
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
        hillMerge(chimeric_hills);

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

        return shrink(new_begin, new_end);
    }

    return true;
}

void Pile::find_repetitive_regions(uint16_t dataset_median) {

    if (!corrected_data_.empty()) {
        corrected_data_.swap(data_);
    }

    if (median_ > 1.42 * dataset_median) {
        dataset_median = std::max(dataset_median, p10_);
    }

    auto slope_regions = find_slopes(1.3);
    if (slope_regions.empty()) {
        if (!corrected_data_.empty()) {
            corrected_data_.swap(data_);
        }
        return;
    }

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
            auto hill_type = hillType(slope_regions[i], slope_regions[j],
                dataset_median, data_, begin_, end_);

            if (hill_type == HillType::kRepeat) {
                uint32_t begin = slope_regions[i].second - 0.336 *
                    (slope_regions[i].second - (slope_regions[i].first >> 1));
                uint32_t end = (slope_regions[j].first >> 1) + 0.336 *
                    (slope_regions[j].second - (slope_regions[j].first >> 1));

                repeat_hills.emplace_back(begin, end);
            }
        }
    }

    // store repetitive hills
    hillMerge(repeat_hills);
    for (const auto& it: repeat_hills) {
        uint32_t hill_begin = std::max(begin_, it.first),
            hill_end = std::min(end_, it.second);
        if (hill_begin < hill_end) {
            hills_.emplace_back(hill_begin, hill_end);
        }
    }

    if (!corrected_data_.empty()) {
        corrected_data_.swap(data_);
    }
}

void Pile::add_repetitive_region(uint32_t begin, uint32_t end) {

    if (begin > data_.size() || end > data_.size()) {
        fprintf(stderr, "[rala::Pile::add_repetitive_region] error: "
            "[begin,end] out of bounds!\n");
        exit(1);
    }

    hills_.emplace_back(begin, end);
}

void Pile::resolve_peculiar_regions() {

    hillMerge(peculiars_);

    auto interval_overlap = [](uint32_t xb, uint32_t xe, uint32_t yb, uint32_t ye) -> uint32_t {
        if (ye < xb || xe < yb) {
            return 0;
        }
        if (xb > yb) {
            return ye - xb;
        }
        return xe - yb;
    };

    std::vector<std::pair<uint32_t, uint32_t>> tmp;
    for (const auto& it: peculiars_) {
        uint32_t l = interval_overlap(it.first, it.second, begin_, end_);
        if (l == 0) {
            continue;
        }

        uint32_t middle = (it.first + it.second) / 2;
        if (l > 0.9 * (it.second - it.first)) {
            tmp.emplace_back(it.first, middle);
            tmp.emplace_back(middle, it.second);
        } else {
            uint32_t ll = interval_overlap(it.first, middle, begin_, end_);
            if (ll == 0) {
                continue;
            }

            uint32_t rl = interval_overlap(middle, it.second, begin_, end_);
            if (rl == 0) {
                continue;
            }

            if (ll > rl) {
                tmp.emplace_back(middle, it.second);
            } else {
                tmp.emplace_back(it.first, middle);
            }
        }
    }
    tmp.swap(peculiars_);
}

void Pile::add_peculiar_region(uint32_t begin, uint32_t end) {

    if (begin > data_.size() || end > data_.size()) {
        fprintf(stderr, "[rala::Pile::add_peculiar_region] error: "
            "[begin,end] out of bounds!\n");
        exit(1);
    }

    peculiars_.emplace_back(begin, end);
}

bool Pile::is_valid_overlap(uint32_t begin, uint32_t end) const {

    begin += begin_;
    end += begin_;
    uint32_t fuzz = 0.042 * (end_ - begin_);

    auto check_hills = [&](const std::vector<std::pair<uint32_t, uint32_t>>& hills) -> bool {
        for (const auto& it: hills) {
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
    };

    return check_hills(hills_) & check_hills(peculiars_);
}

std::string Pile::to_json() const {

    std::stringstream ss;
    ss << "\"" << id_ << "\":{";

    ss << "\"y\":[";
    for (uint32_t i = 0; i < data_.size(); ++i) {
        ss << data_[i];
        if (i < data_.size() - 1) {
            ss << ",";
        }
    }
    ss << "],";

    ss << "\"~y\":[";
    for (uint32_t i = 0; i < corrected_data_.size(); ++i) {
        ss << corrected_data_[i];
        if (i < corrected_data_.size() - 1) {
            ss << ",";
        }
    }
    ss << "],";

    ss << "\"b\":" << begin_ << ",";
    ss << "\"e\":" << end_ << ",";

    ss << "\"h\":[";
    for (uint32_t i = 0; i < hills_.size(); ++i) {
        ss << hills_[i].first << "," << hills_[i].second;
        if (i < hills_.size() - 1) {
            ss << ",";
        }
    }

    if (!hills_.empty() && !peculiars_.empty()) {
        ss << ",";
    }

    for (uint32_t i = 0; i < peculiars_.size(); ++i) {
        ss << peculiars_[i].first << "," << peculiars_[i].second;
        if (i < peculiars_.size() - 1) {
            ss << ",";
        }
    }
    ss << "],";

    ss << "\"m\":" << median_ << ",";
    ss << "\"p10\":" << p10_;
    ss << "}";

    return ss.str();
}

}
