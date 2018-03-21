/*!
 * @file sequence.cpp
 *
 * @brief Sequence class source file
 */

#include <stdlib.h>

#include "sequence.hpp"

namespace rala {

std::unique_ptr<Sequence> createSequence(const std::string& name,
    const std::string& data) {

    if (name.empty()) {
        fprintf(stderr, "[rala::createSequence] error: empty name!\n");
        exit(1);
    }
    if (data.empty()) {
        fprintf(stderr, "[rala::createSequence] error: empty data!\n");
        exit(1);
    }

    return std::unique_ptr<Sequence>(new Sequence(name.c_str(), name.size(),
        data.c_str(), data.size()));
}

Sequence::Sequence(const char* name, uint32_t name_length, const char* data,
    uint32_t data_length)
        : name_(name, name_length), data_(data, data_length), reverse_complement_() {
}

Sequence::Sequence(const char* name, uint32_t name_length, const char* data,
    uint32_t data_length, const char*, uint32_t)
        : Sequence(name, name_length, data, data_length) {
}

void Sequence::trim(uint32_t begin, uint32_t end) {
    data_ = data_.substr(begin, end - begin);
    if (!reverse_complement_.empty()) {
        create_reverse_complement();
    }
}

void Sequence::create_reverse_complement() {

    reverse_complement_.clear();

    for (int32_t i = data_.size() - 1; i >= 0; --i) {
        switch (data_[i]) {
            case 'A':
                reverse_complement_ += 'T';
                break;
            case 'T':
                reverse_complement_ += 'A';
                break;
            case 'C':
                reverse_complement_ += 'G';
                break;
            case 'G':
                reverse_complement_ += 'C';
                break;
            default:
                reverse_complement_ += data_[i];
                break;
        }
    }
}

}
