/*!
 * @file overlap.hpp
 *
 * @brief Overlap class header file
 */

#include "overlap.hpp"

namespace RALAY {

std::unique_ptr<Overlap> createOverlap(uint32_t id, const double* values,
    uint32_t values_length) {

    return std::unique_ptr<Overlap>(new Overlap(id, values, values_length));
}

Overlap::Overlap(uint32_t id, const double* values, uint32_t values_length) :
        id_(id), error_(values[2]), minmers_(values[3]), a_id_(values[0]),
        a_rc_(values[4]), a_begin_(values[5]), a_end_(values[6]),
        a_length_(values[7]), b_id_(values[1]), b_rc_(values[8]),
        b_begin_(values[9]), b_end_(values[10]), b_length_(values[11]) {
}

Overlap::~Overlap() {
}

}
