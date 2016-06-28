/*!
 * @file read.cpp
 *
 * @brief Read class source file
 */

#include <assert.h>

#include "read.hpp"

namespace RALAY {

std::unique_ptr<Read> createRead(uint32_t id, const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length, const char* quality,
    uint32_t quality_length) {

    assert(sequence_length > 0);
    assert(sequence_length == quality_length);

    return std::unique_ptr<Read>(new Read(id, name, name_length, sequence,
        sequence_length, quality, quality_length));
}

Read::Read(uint32_t id, const char* name, uint32_t name_length, const char* sequence,
    uint32_t sequence_length, const char* quality, uint32_t quality_length) :
        id_(id), name_(name, name_length), sequence_(sequence, sequence_length),
        quality_(quality, quality_length) {
}

}
