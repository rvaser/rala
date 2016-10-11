/*!
 * @file read.hpp
 *
 * @brief Read class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>

#include "bioparser/src/bioparser.hpp"

namespace RALA {

class Read;
std::unique_ptr<Read> createRead(uint32_t id, const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length, const char* quality,
    uint32_t quality_length);

class Read {
    public:
        ~Read() {};

        uint32_t id() const {
            return id_;
        }

        const std::string& name() const {
            return name_;
        }

        const std::string& sequence() const {
            return sequence_;
        }

        const std::string& quality() const {
            return quality_;
        }

        const std::string& rc() {
            if (rc_.size() != sequence_.size()) create_rc();
            return rc_;
        }

        void trim_sequence(uint32_t begin, uint32_t end);

        friend std::unique_ptr<Read> createRead(uint32_t id, const char* name,
            uint32_t name_length, const char* sequence, uint32_t sequence_length,
            const char* quality, uint32_t quality_length);

        friend BIOPARSER::FastqReader<Read>;

    private:
        Read(uint32_t id, const char* name, uint32_t name_length, const char* sequence,
            uint32_t sequence_length, const char* quality, uint32_t quality_length);
        Read(const Read&) = delete;
        const Read& operator=(const Read&) = delete;

        void create_rc();

        uint32_t id_;
        std::string name_;
        std::string sequence_;
        std::string quality_;
        std::string rc_;
};

}
