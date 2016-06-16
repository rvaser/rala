/*!
 * @file Read.hpp
 *
 * @brief Read class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>

namespace RALAY {

class Read;
std::unique_ptr<Read> createRead(uint32_t id, const std::string& sequence,
    const std::string& quality);

class Read {
public:

    ~Read() {};

    int id() const {
        return id_;
    }

    const std::string& sequence() const {
        return sequence_;
    }

    const std::string& quality() const {
        return quality_;
    }

    friend std::unique_ptr<Read> createRead(uint32_t id, const std::string& sequence,
        const std::string& quality);

private:

    Read(uint32_t id, const std::string& sequence, const std::string& quality);
    Read(const Read&) = delete;
    const Read& operator=(const Read&) = delete;

    uint32_t id_;
    std::string sequence_;
    std::string quality_;
};

}
