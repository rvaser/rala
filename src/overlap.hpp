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
std::unique_ptr<Overlap> createOverlap(uint32_t id, uint32_t a_id, uint32_t b_id,
    double error, uint32_t minmers, uint32_t a_rc, uint32_t a_begin, uint32_t a_end,
    uint32_t a_length, uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length);


class Overlap {
public:

    ~Overlap();

    uint32_t id() const {
        return id_;
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

    uint32_t b_rc() const {
        return b_rc_;
    }

    uint32_t b_id() const {
        return b_id_;
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

    double quality() const {
        return quality_;
    }

    uint32_t length() const {
        return length_;
    }

    uint32_t matching_bases() const {
        return matching_bases_;
    }

    // returns whether the new overlap is valid
    bool update(uint32_t a_begin, uint32_t a_end, uint32_t b_begin, uint32_t b_end);

    friend std::unique_ptr<Overlap> createOverlap(uint32_t id, uint32_t a_id, uint32_t b_id,
        double error, uint32_t minmers, uint32_t a_rc, uint32_t a_begin, uint32_t a_end,
        uint32_t a_length, uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length);

    friend BIOPARSER::MhapReader<Overlap>;

private:

    Overlap(uint32_t id, uint32_t a_id, uint32_t b_id, double error, uint32_t minmers,
        uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
        uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length);
    Overlap(const Overlap&) = delete;
    const Overlap& operator=(const Overlap&) = delete;

    uint32_t id_;

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

    double quality_;
    uint32_t length_;
    uint32_t matching_bases_;
};

}
