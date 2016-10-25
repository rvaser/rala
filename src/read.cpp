/*!
 * @file read.cpp
 *
 * @brief Read class source file
 */

#include <assert.h>

#include "read.hpp"

namespace RALA {

std::unique_ptr<Read> createRead(uint32_t id, const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length, const char* quality,
    uint32_t quality_length) {

    assert(sequence_length > 0);
    assert(sequence_length == quality_length);

    return std::unique_ptr<Read>(new Read(id, name, name_length, sequence,
        sequence_length, quality, quality_length));
}

Read::Read(uint32_t id, const char* name, uint32_t name_length, const char* sequence,
    uint32_t sequence_length, const char* quality, uint32_t quality_length)
        : id_(id), name_(name, name_length), sequence_(sequence, sequence_length),
        quality_(quality, quality_length), rc_() {
}

void Read::trim_sequence(uint32_t begin, uint32_t end) {
    sequence_ = sequence_.substr(begin, end - begin);
    quality_ = quality_.substr(begin, end - begin);
}

void Read::create_rc() {

    rc_.clear();
    for (int32_t i = sequence_.size() - 1; i >= 0; --i) {
        char c = sequence_[i];
        switch (c) {
            case 'A':
                c = 'T';
                break;
            case 'T':
                c = 'A';
                break;
            case 'C':
                c = 'G';
                break;
            case 'G':
                c = 'C';
                break;
            default:
                break;
        }
        rc_ += c;
    }
}

}
