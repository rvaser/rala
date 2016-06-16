/*!
 * @file Read.cpp
 *
 * @brief Read class source file
 */

#include <assert.h>

#include "read.hpp"

namespace RALAY {

std::unique_ptr<Read> createRead(uint32_t id, const std::string& sequence,
    const std::string& quality) {

    assert(sequence.size() != 0);
    assert(sequence.size() == quality.size());

    return std::unique_ptr<Read>(new Read(id, sequence, quality));
}

Read::Read(uint32_t id, const std::string& sequence, const std::string& quality) :
        id_(id), sequence_(sequence), quality_(quality) {
}

}
