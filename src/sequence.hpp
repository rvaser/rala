/*!
 * @file sequence.hpp
 *
 * @brief Sequence class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>

namespace bioparser {
    template<class T>
    class FastaParser;

    template<class T>
    class FastqParser;
}

namespace rala {

class Sequence;
std::unique_ptr<Sequence> createSequence(const std::string& name,
    const std::string& data);

class Sequence {
public:
    ~Sequence() {};

    std::uint64_t id() const {
        return id_;
    }

    const std::string& name() const {
        return name_;
    }

    const std::string& data() const {
        return data_;
    }

    const std::string& reverse_complement() {
        if (reverse_complement_.size() != data_.size()) {
            create_reverse_complement();
        }
        return reverse_complement_;
    }

    void trim(uint32_t begin, uint32_t end);

    friend std::unique_ptr<Sequence> createSequence(const std::string& name,
        const std::string& data);

    friend bioparser::FastaParser<Sequence>;
    friend bioparser::FastqParser<Sequence>;

    static std::uint64_t num_objects;
private:
    Sequence(const char* name, uint32_t name_length, const char* data,
        uint32_t data_length);
    Sequence(const char* name, uint32_t name_length, const char* data,
        uint32_t data_length, const char* quality, uint32_t quality_length);
    Sequence(const Sequence&) = delete;
    const Sequence& operator=(const Sequence&) = delete;

    void create_reverse_complement();

    std::uint64_t id_;
    std::string name_;
    std::string data_;
    std::string reverse_complement_;
};

}
