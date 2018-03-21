/*!
 * @file pile.hpp
 *
 * @brief Pile class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>

namespace rala {

class Overlap;

class Pile;
std::unique_ptr<Pile> createPile(uint64_t id, uint32_t sequence_length);

class Pile {
public:
    ~Pile() {};

    uint64_t id() const {
        return id_;
    }

    /*!
     * @brief Returns begin_ of the valid interval [begin_, end_>
     */
    uint32_t begin() const {
        return begin_;
    }

    /*!
     * @brief Returns end_ of the valid interval [begin_, end_>
     */
    uint32_t end() const {
        return end_;
    }

    uint16_t median() const {
        return median_;
    };

    void find_median();

    const std::vector<uint16_t>& data() const {
        if (corrected_data_.empty()) {
            return data_;
        }
        return corrected_data_;
    }

    /*!
     * @brief Adds overlaps to data_
     */
    void add_layers(std::vector<uint32_t>& overlap_bounds);

    /*!
     * @brief Sets values of data_ outside the interval [begin, end> to zeroes
     * and updates begin_, end_ accordingly
     */
    bool shrink(uint32_t begin, uint32_t end);

    /*!
     * @brief Corrects data_ with other piles which have overlapping regions
     */
    void correct(const std::vector<std::shared_ptr<Overlap>>& overlaps,
        const std::vector<std::unique_ptr<Pile>>& piles);

    /*!
     * @brief Locates region in data_ with values greater or equal to predefined
     * coverage; updates begin_, end_ and data_ accordingly;
     * if there is no such region (with valid coverage and longer than 1000),
     * false is returned
     */
    bool find_valid_region();

    /*!
     * @brief Locates chimeric regions in data_;
     * if any such region is found, both begin_ and end_ are set to the longest
     * continuous region of the sequence and data_ is updated accordingly;
     * if the new area is shorter than 1000, false is returned
     */
    bool find_chimeric_regions(uint16_t dataset_median);

    /*!
     * @brief Locates regions in data_ which ought to be repetitive in the
     * genome and stores them in hills_
     */
    void find_repetitive_regions(uint16_t dataset_median);

    /*!
     * @brief Checks whether overlap [begin, end> is valid with respect to
     * hills_ which indicate repetitive regions of the genome
     */
    bool is_valid_overlap(uint32_t begin, uint32_t end) const;

    /*!
     * @brief Prints data_ in csv format to path
     */
    void print_csv(std::string path, uint16_t dataset_median = 0) const;

    friend std::unique_ptr<Pile> createPile(uint64_t id, uint32_t sequence_length);
private:
    Pile(uint64_t id, uint32_t sequence_length);
    Pile(const Pile&) = delete;
    const Pile& operator=(const Pile&) = delete;

    std::vector<std::pair<uint32_t, uint32_t>> find_slopes(double q);

    uint64_t id_;
    uint32_t begin_;
    uint32_t end_;
    uint16_t p10_;
    uint16_t median_;
    std::vector<uint16_t> data_;
    std::vector<uint16_t> corrected_data_;
    std::vector<std::pair<uint32_t, uint32_t>> hills_;
};

}
