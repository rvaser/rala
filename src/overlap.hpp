/*!
 * @file overlap.cpp
 *
 * @brief Overlap class source file
 */

#pragma once

#include <stdint.h>
#include <memory>

namespace bioparser {
    template<class T>
    class MhapReader;

    template<class T>
    class PafReader;
}

namespace rala {

enum class OverlapType {
    kX, // internal
    kA, // b contained
    kB, // a contained
    kAB, // suffix prefix
    kBA // prefix suffix
};

class Overlap {
public:
    ~Overlap();

    uint64_t id() const {
        return id_;
    }

    uint32_t a_id() const {
        return a_id_;
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

    uint32_t b_begin() const {
        return b_begin_;
    }

    uint32_t b_end() const {
        return b_end_;
    }

    uint32_t b_length() const {
        return b_length_;
    }

    uint32_t length() const {
        return length_;
    }

    uint32_t orientation() const {
        return orientation_;
    }

    OverlapType type() const;

    bool update(uint32_t a_trimmed_begin, uint32_t a_trimmed_end,
        uint32_t b_trimmed_begin, uint32_t b_trimmed_end);

    friend bioparser::MhapReader<Overlap>;
    friend bioparser::PafReader<Overlap>;
private:
    Overlap(uint64_t id, uint32_t a_id, uint32_t b_id, double error,
        uint32_t minmers, uint32_t a_rc, uint32_t a_begin, uint32_t a_end,
        uint32_t a_length, uint32_t b_rc, uint32_t b_begin, uint32_t b_end,
        uint32_t b_length);
    Overlap(uint64_t id, const char* a_name, uint32_t a_name_length,
        uint32_t a_length, uint32_t a_begin, uint32_t a_end,
        char orientation, const char* b_name, uint32_t b_name_length,
        uint32_t b_length, uint32_t b_begin, uint32_t b_end,
        uint32_t matching_bases, uint32_t overlap_length, uint32_t quality);
    Overlap(const Overlap&) = delete;
    const Overlap& operator=(const Overlap&) = delete;

    uint64_t id_;
    uint32_t a_id_;
    uint32_t a_begin_;
    uint32_t a_end_;
    uint32_t a_length_;
    uint32_t b_id_;
    uint32_t b_begin_;
    uint32_t b_end_;
    uint32_t b_length_;
    uint32_t length_;
    uint32_t orientation_;
};

}
