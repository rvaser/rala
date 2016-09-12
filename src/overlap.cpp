/*!
 * @file overlap.hpp
 *
 * @brief Overlap class header file
 */

#include "overlap.hpp"

namespace RALA {

std::unique_ptr<Overlap> createOverlap(uint32_t id, uint32_t a_id, uint32_t b_id,
    double error, uint32_t minmers, uint32_t a_rc, uint32_t a_begin, uint32_t a_end,
    uint32_t a_length, uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length) {

    return std::unique_ptr<Overlap>(new Overlap(id, a_id, b_id, error, minmers,
        a_rc, a_begin, a_end, a_length, b_rc, b_begin, b_end, b_length));
}

Overlap::Overlap(uint32_t id, uint32_t a_id, uint32_t b_id, double error, uint32_t minmers,
    uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
    uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length) :
        id_(id), a_id_(a_id - 1), a_rc_(a_rc), a_begin_(a_begin), a_end_(a_end), a_length_(a_length),
        b_id_(b_id - 1), b_rc_(b_rc), b_begin_(b_begin), b_end_(b_end), b_length_(b_length),
        quality_(error), length_(std::max(a_end - a_begin, b_end - b_begin)),
        matching_bases_((uint32_t) (length_ * quality_ + 0.499)), type_() {
}

Overlap::~Overlap() {
}

bool Overlap::update(uint32_t a_trimmed_begin, uint32_t a_trimmed_end,
    uint32_t b_trimmed_begin, uint32_t b_trimmed_end) {

    assert(a_trimmed_begin <= a_length_);
    assert(a_trimmed_end <= a_length_);
    assert(b_trimmed_begin <= b_length_);
    assert(b_trimmed_end <= b_length_);

    if (a_rc_) {
        uint32_t tmp = a_trimmed_begin;
        a_trimmed_begin = a_length_ - a_trimmed_end;
        a_trimmed_end = a_length_ - tmp;
    }
    if (b_rc_) {
        uint32_t tmp = b_trimmed_begin;
        b_trimmed_begin = b_length_ - b_trimmed_end;
        b_trimmed_end = b_length_ - tmp;
    }

    uint32_t a_new_begin = a_begin_ + (b_begin_ < b_trimmed_begin ? b_trimmed_begin - b_begin_ : 0);
    uint32_t a_new_end = a_end_ - (b_end_ > b_trimmed_end ? b_end_ - b_trimmed_end : 0);
    uint32_t b_new_begin = b_begin_ + (a_begin_ < a_trimmed_begin ? a_trimmed_begin - a_begin_ : 0);
    uint32_t b_new_end = b_end_ - (a_end_ > a_trimmed_end ? a_end_ - a_trimmed_end : 0);

    a_new_begin = (a_new_begin < a_trimmed_begin ? a_trimmed_begin : a_new_begin) - a_trimmed_begin;
    a_new_end = (a_new_end > a_trimmed_end ? a_trimmed_end : a_new_end) - a_trimmed_begin;
    b_new_begin = (b_new_begin < b_trimmed_begin ? b_trimmed_begin : b_new_begin) - b_trimmed_begin;
    b_new_end = (b_new_end > b_trimmed_end ? b_trimmed_end : b_new_end) - b_trimmed_begin;

    if (a_new_begin >= a_new_end || b_new_begin >= b_new_end) {
        return false;
    }

    double ratio = ((a_new_end - a_new_begin) + (b_new_end - b_new_begin)) /
        (double) ((a_end_ - a_begin_) + (b_end_ - b_begin_));

    length_ = (uint32_t) (length_ * ratio + 0.499);
    matching_bases_ = (uint32_t) (matching_bases_ * ratio + 0.499);
    quality_ = (matching_bases_ / (double) length_);

    a_begin_ = a_new_begin;
    a_end_ = a_new_end;
    a_length_ = a_trimmed_end - a_trimmed_begin;

    b_begin_ = b_new_begin;
    b_end_ = b_new_end;
    b_length_ = b_trimmed_end - b_trimmed_begin;

    return true;
}

}
