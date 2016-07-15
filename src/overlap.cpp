/*!
 * @file overlap.hpp
 *
 * @brief Overlap class header file
 */

#include "overlap.hpp"

namespace RALAY {

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
        matching_bases_((uint32_t) (length_ * quality_ + 0.499)) {
}

Overlap::~Overlap() {
}

}
