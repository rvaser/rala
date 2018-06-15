/*!
 * @file overlap.cpp
 *
 * @brief Overlap class source file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <unordered_map>

namespace bioparser {
    template<class T>
    class MhapParser;

    template<class T>
    class PafParser;
}

namespace rala {

class Pile;
class Graph;

enum class OverlapType {
    kX, // bad overlap
    kA, // b contained
    kB, // a contained
    kAB, // suffix prefix
    kBA // prefix suffix
};

class Overlap {
public:
    ~Overlap();

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

    bool transmute(const std::vector<std::unique_ptr<Pile>>& piles,
        const std::unordered_map<std::string, uint64_t>& name_to_id);

    bool trim(const std::vector<std::unique_ptr<Pile>>& piles);

    OverlapType type() const;

    friend bioparser::MhapParser<Overlap>;
    friend bioparser::PafParser<Overlap>;
    friend Graph;
private:
    Overlap(uint64_t a_id, uint64_t b_id, double error, uint32_t minmers,
        uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
        uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length);
    Overlap(const char* a_name, uint32_t a_name_length,
        uint32_t a_length, uint32_t a_begin, uint32_t a_end,
        char orientation, const char* b_name, uint32_t b_name_length,
        uint32_t b_length, uint32_t b_begin, uint32_t b_end,
        uint32_t matching_bases, uint32_t overlap_length, uint32_t quality);
    Overlap(const Overlap&) = delete;
    const Overlap& operator=(const Overlap&) = delete;

    std::string a_name_;
    uint64_t a_id_;
    uint32_t a_begin_;
    uint32_t a_end_;
    uint32_t a_length_;
    std::string b_name_;
    uint64_t b_id_;
    uint32_t b_begin_;
    uint32_t b_end_;
    uint32_t b_length_;
    uint32_t length_;
    uint32_t orientation_;
    bool is_transmuted_;
};

}
