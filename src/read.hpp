/*!
 * @file read.hpp
 *
 * @brief Read class header file
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <string>
#include <deque>

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
std::unique_ptr<ReadInfo> copyReadInfo(std::shared_ptr<ReadInfo> other);

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

    /*!
     * @brief Calculates coverage median, call before coverage_median()
     */
    void find_coverage_median();

    const std::vector<uint16_t>& coverage_graph() const {
        return coverage_graph_;
    }

    /*!
     * @brief Adds overlaps to coverage_graph_
     */
    void update_coverage_graph(std::vector<uint32_t>& shrunken_overlaps);

    /*!
     * @brief Reduces coverage_graph_ by setting values outside the interval
     * [begin, end> to zeroes
     */
    bool reduce_coverage_graph(uint32_t begin, uint32_t end);

    /*!
     * @biref Smooths coverage graph with 1D average filter
     */
    void smooth_coverage_graph();

    /*!
     * @brief Correct coverage graph with other graph over overlapping region
     */
    void correct_coverage_graph(uint32_t region_begin, uint32_t region_end,
        std::shared_ptr<ReadInfo> other, uint32_t other_region_begin,
        uint32_t other_region_end, uint32_t orientation);

    /*!
     * @brief Locates region in coverage_graph_ with values greater or equal to
     * predefined value;
     * updates begin_, end_ and coverage_graph_ accordingly;
     * if there is no such region (with valid coverage and longer than 1000),
     * object is invalidated
     */
    bool find_valid_region();

    /*!
     * @brief Locates pits in coverage_graph_ which ought to indicate that the
     * read is chimeric;
     * if a pit is found, both begin_ and end_ are set to the longest continuous
     * region of the read and coverage_graph_ is updated accordingly;
     * if the new area is shorter than 500, object is invalidated;
     * if no pits are found, false is returned
     */
    bool find_coverage_pits(uint16_t dataset_median);

    /*!
     * @brief Returns regions of read which ought to be repetitive in the genome
     */
    const std::vector<std::pair<uint32_t, uint32_t>>& coverage_hills() const {
        return coverage_hills_;
    }

    /*!
     * @brief Locates regions in coverage_graph_ which ought to be repetitive
     * in the genome and stores them in coverage_hills_
     */
    void find_coverage_hills(uint16_t dataset_median);

    /*!
     * @brief Print coverage_graph_ in csv format to path
     */
    void print_csv(std::string path, uint16_t dataset_median = 0) const;

    friend std::unique_ptr<ReadInfo> createReadInfo(uint64_t id,
        uint32_t read_length);

    friend std::unique_ptr<ReadInfo> copyReadInfo(
        std::shared_ptr<ReadInfo> read_info);
private:
    ReadInfo(uint64_t id, uint32_t read_length);
    ReadInfo(const ReadInfo&) = default;
    const ReadInfo& operator=(const ReadInfo&) = delete;

    static void coverage_window_add(
        std::deque<std::pair<int32_t, int32_t>>& window, int32_t value,
        int32_t position);
    static void coverage_window_update(
        std::deque<std::pair<int32_t, int32_t>>& window, int32_t position);

    uint64_t id_;
    uint32_t begin_;
    uint32_t end_;
    uint16_t coverage_median_;
    std::vector<uint16_t> coverage_graph_;
    std::vector<std::pair<uint32_t, uint32_t>> coverage_hills_;
};

}
