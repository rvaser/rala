/*!
 * @file overlap.cpp
 *
 * @brief Overlap class source file
 */

#pragma once

#include <stdint.h>
#include <memory>

#include "bioparser/src/bioparser.hpp"

namespace RALAY {

class Overlap;
std::unique_ptr<Overlap> createOverlap(uint32_t id, const double* values,
    uint32_t values_length);


class Overlap {
public:

    ~Overlap();

    uint32_t id() const {
        return id_;
    }

    double error() const {
        return error_;
    }

    uint32_t minmers() const {
        return minmers_;
    }

    uint32_t a_id() const {
        return a_id_;
    }

    uint32_t a_rc() const {
        return a_rc_;
    }

    uint32_t a_begin() const {
        return a_begin_;
    }

    uint32_t a_end() const {
        return a_end_;
    }

    uint32_t a_length() const {
        return a_length_;
    }

    uint32_t b_id() const {
        return b_id_;
    }

    uint32_t b_rc() const {
        return b_rc_;
    }

    uint32_t b_begin() const {
        return b_begin_;
    }

    uint32_t b_end() const {
        return b_end_;
    }

    uint32_t b_length() const {
        return b_length_;
    }

    friend std::unique_ptr<Overlap> createOverlap(uint32_t id, const double* values,
        uint32_t values_length);

    friend BIOPARSER::MhapReader<Overlap>;

private:

    Overlap(uint32_t id, const double* values, uint32_t values_length);
    Overlap(const Overlap&) = delete;
    const Overlap& operator=(const Overlap&) = delete;

    uint32_t id_;
    double error_;
    uint32_t minmers_;

    uint32_t a_id_;
    uint32_t a_rc_;
    uint32_t a_begin_;
    uint32_t a_end_;
    uint32_t a_length_;

    uint32_t b_id_;
    uint32_t b_rc_;
    uint32_t b_begin_;
    uint32_t b_end_;
    uint32_t b_length_;
};

}
