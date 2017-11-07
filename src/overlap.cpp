/*!
 * @file overlap.hpp
 *
 * @brief Overlap class header file
 */

#include <assert.h>

#include "overlap.hpp"

namespace rala {

Overlap::Overlap(uint64_t id, uint32_t a_id, uint32_t b_id, double error,
    uint32_t minmers, uint32_t a_rc, uint32_t a_begin, uint32_t a_end,
    uint32_t a_length, uint32_t b_rc, uint32_t b_begin, uint32_t b_end,
    uint32_t b_length)
        : id_(id), a_id_(a_id - 1), a_begin_(a_begin), a_end_(a_end),
        a_length_(a_length), b_id_(b_id - 1), b_begin_(b_begin), b_end_(b_end),
        b_length_(b_length), length_(std::max(a_end - a_begin, b_end - b_begin)),
        orientation_(a_rc == b_rc ? 0 : 1) {
}

Overlap::Overlap(uint64_t id, const char* a_name, uint32_t a_name_length,
    uint32_t a_length, uint32_t a_begin, uint32_t a_end, char orientation,
    const char* b_name, uint32_t b_name_length, uint32_t b_length,
    uint32_t b_begin, uint32_t b_end, uint32_t matching_bases,
    uint32_t overlap_length, uint32_t quality)
        : id_(id), a_id_(atoi(a_name) - 1), a_begin_(a_begin), a_end_(a_end),
        a_length_(a_length), b_id_(atoi(b_name) - 1), b_begin_(b_begin),
        b_end_(b_end), b_length_(b_length), length_(overlap_length),
        orientation_(orientation == '+' ? 0 : 1) {
}

Overlap::~Overlap() {
}

bool Overlap::update(uint32_t a_trimmed_begin, uint32_t a_trimmed_end,
    uint32_t b_trimmed_begin, uint32_t b_trimmed_end) {

    assert(a_trimmed_begin <= a_length_);
    assert(a_trimmed_end <= a_length_);
    assert(b_trimmed_begin <= b_length_);
    assert(b_trimmed_end <= b_length_);

    if (a_begin_ >= a_trimmed_end || a_end_ <= a_trimmed_begin) {
        return false;
    }
    if (b_begin_ >= b_trimmed_end || b_end_ <= b_trimmed_begin) {
        return false;
    }

    uint32_t a_new_begin = a_begin_ + (b_begin_ < b_trimmed_begin ?
        b_trimmed_begin - b_begin_ : 0);
    uint32_t a_new_end = a_end_ - (b_end_ > b_trimmed_end ?
        b_end_ - b_trimmed_end : 0);
    if (a_new_begin >= a_trimmed_end || a_new_end <= a_trimmed_begin) {
        return false;
    }

    uint32_t b_new_begin = b_begin_ + (a_begin_ < a_trimmed_begin ?
        a_trimmed_begin - a_begin_ : 0);
    uint32_t b_new_end = b_end_ - (a_end_ > a_trimmed_end ?
        a_end_ - a_trimmed_end : 0);
    if (b_new_begin >= b_trimmed_end || b_new_end <= b_trimmed_begin) {
        return false;
    }

    a_new_begin = (a_new_begin < a_trimmed_begin ? a_trimmed_begin :
        a_new_begin) - a_trimmed_begin;
    a_new_end = (a_new_end > a_trimmed_end ? a_trimmed_end : a_new_end) -
        a_trimmed_begin;
    b_new_begin = (b_new_begin < b_trimmed_begin ? b_trimmed_begin :
        b_new_begin) - b_trimmed_begin;
    b_new_end = (b_new_end > b_trimmed_end ? b_trimmed_end : b_new_end) -
        b_trimmed_begin;

    if (a_new_begin >= a_new_end || b_new_begin >= b_new_end) {
        return false;
    }

    a_begin_ = a_new_begin;
    a_end_ = a_new_end;
    a_length_ = a_trimmed_end - a_trimmed_begin;

    b_begin_ = b_new_begin;
    b_end_ = b_new_end;
    b_length_ = b_trimmed_end - b_trimmed_begin;

    length_ = std::max(a_end_ - a_begin_, b_end_ - b_begin_);

    return true;
}

OverlapType Overlap::type() const {

    uint32_t a_begin = a_begin_;
    uint32_t a_end = a_end_;
    uint32_t b_begin = orientation_ == 0 ? b_begin_ : b_length_ - b_end_;
    uint32_t b_end = orientation_ == 0 ? b_end_ : b_length_ - b_begin_;

    uint32_t overhang = std::min(a_begin, b_begin) + std::min(a_length_ -
        a_end, b_length_ - b_end);

    if (a_end - a_begin < (a_end - a_begin + overhang) * 0.875 ||
        b_end - b_begin < (b_end - b_begin + overhang) * 0.875) {
        return OverlapType::kX;
    }
    if (a_begin <= b_begin && (a_length_ - a_end) <= (b_length_ - b_end)) {
        return OverlapType::kB;
    }
    if (a_begin >= b_begin && (a_length_ - a_end) >= (b_length_ - b_end)) {
        return OverlapType::kA;
    }

    auto absolute_difference = [](uint32_t a, uint32_t b) -> uint32_t {
        return a > b ? (a - b) : (b - a);
    };

    uint32_t overlap_length = std::max(a_end_ - a_begin_, b_end_ - b_begin_);

    if (absolute_difference(a_end_ - a_begin_, b_end_ - b_begin_) < overlap_length * 0.01) {
        uint32_t min_extension = 0.05 * std::max(a_length_, b_length_);

        if (absolute_difference(a_begin, b_begin) < min_extension) {
            if ((a_length_ - a_end) >= (b_length_ - b_end)) {
                return OverlapType::kA;
            } else {
                return OverlapType::kB;
            }
        }
        if (absolute_difference((a_length_ - a_end), (b_length_ - b_end)) < min_extension) {
            if (a_begin >= b_begin) {
                return OverlapType::kA;
            } else {
                return OverlapType::kB;
            }
        }
    }

    if (a_begin > b_begin) {
        return OverlapType::kAB;
    }
    return OverlapType::kBA;
}

}
