/*!
 * @file read.hpp
 *
 * @brief Read class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>

namespace bioparser {
    template<class T>
    class FastaReader;

    template<class T>
    class FastqReader;
}

namespace rala {

class Read {
public:
    ~Read() {};

    uint64_t id() const {
        return id_;
    }

    const std::string& name() const {
        return name_;
    }

    const std::string& sequence() const {
        return sequence_;
    }

    const std::string& reverse_complement() {
        if (reverse_complement_.size() != sequence_.size()) {
            create_reverse_complement();
        }
        return reverse_complement_;
    }

    void update(uint32_t begin, uint32_t end);

    friend bioparser::FastaReader<Read>;
    friend bioparser::FastqReader<Read>;
private:
    Read(uint64_t id, const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length);
    Read(uint64_t id, const char* name, uint32_t name_length,
        const char* sequence, uint32_t sequence_length,
        const char* quality, uint32_t quality_length);
    Read(const Read&) = delete;
    const Read& operator=(const Read&) = delete;

    void create_reverse_complement();

    uint64_t id_;
    std::string name_;
    std::string sequence_;
    std::string reverse_complement_;
};

class ReadInfo;
std::unique_ptr<ReadInfo> createReadInfo(uint64_t id, uint32_t read_length);

class ReadInfo {
public:
    ~ReadInfo() {};

    uint64_t id() const {
        return id_;
    }

    /*!
     * @brief Returns begin_ of the valid coverage interval [begin_, end_>
     */
    uint32_t begin() const {
        return begin_;
    }

    /*!
     * @brief Returns end_ of the valid coverage interval [begin_, end_>
     */
    uint32_t end() const {
        return end_;
    }

    uint16_t coverage_median() const {
        return coverage_median_;
    };

    const std::vector<uint16_t>& coverage_graph() const {
        if (corrected_coverage_graph_.empty()) {
            return coverage_graph_;
        }
        return corrected_coverage_graph_;
    }

    void add_coverage_hill(uint32_t begin, uint32_t end) {
        coverage_hills_.emplace_back(begin, end);
    }

    /*!
     * @brief Adds overlaps to coverage_graph_
     */
    void update_coverage_graph(std::vector<uint32_t>& shrunken_overlaps);

    /*!
     * @brief Sets values of coverage_graph_ outside the interval [begin, end>
     * to zeroes and updates begin_, end_ accordingly
     */
    bool shrink_coverage_graph(uint32_t begin, uint32_t end);

    /*!
     * @brief Correct coverage graph with other graph over overlapping region
     */
    void correct_coverage_graph(const std::vector<uint32_t>& shrunken_overlaps,
        const std::vector<std::unique_ptr<ReadInfo>>& read_infos);

    /*!
     * @brief Calculates coverage median, call before coverage_median()
     */
    void find_coverage_median();

    /*!
     * @brief Locates region in coverage_graph_ with values greater or equal to
     * predefined coverage; updates begin_, end_ and coverage_graph_ accordingly;
     * if there is no such region (with valid coverage and longer than 1000),
     * false is returned
     */
    bool find_valid_region();

    /*!
     * @brief Locates chimeric regions in coverage_graph_;
     * if any such region is found, both begin_ and end_ are set to the longest
     * continuous region of the read and coverage_graph_ is updated accordingly;
     * if the new area is shorter than 1000, false is returned
     */
    bool find_chimeric_region(uint16_t dataset_median);

    /*!
     * @brief Locates regions in coverage_graph_ which ought to be repetitive
     * in the genome and stores them in coverage_hills_
     */
    void find_repetitive_region(uint16_t dataset_median);

    /*!
     * @brief Checks whether overlap [begin, end> is valid with respect to
     * coverage_hills_ which indicate repetitive regions of the genome
     */
    bool is_valid_overlap(uint32_t begin, uint32_t end) const;

    /*!
     * @brief Print coverage_graph_ in csv format to path
     */
    void print_csv(std::string path, uint16_t dataset_median = 0) const;

    friend std::unique_ptr<ReadInfo> createReadInfo(uint64_t id,
        uint32_t read_length);
private:
    ReadInfo(uint64_t id, uint32_t read_length);
    ReadInfo(const ReadInfo&) = delete;
    const ReadInfo& operator=(const ReadInfo&) = delete;

    std::vector<std::pair<uint32_t, uint32_t>> find_coverage_slopes(double q);

    uint64_t id_;
    uint32_t begin_;
    uint32_t end_;
    uint16_t coverage_p10_;
    uint16_t coverage_median_;
    std::vector<uint16_t> coverage_graph_;
    std::vector<uint16_t> corrected_coverage_graph_;
    std::vector<std::pair<uint32_t, uint32_t>> coverage_hills_;
};

}
