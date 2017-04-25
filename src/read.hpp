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

#include "bioparser/src/bioparser.hpp"

namespace rala {

class Read;
std::unique_ptr<Read> createRead(uint64_t id, const char* name, uint32_t name_length,
    const char* sequence, uint32_t sequence_length, const char* quality,
    uint32_t quality_length);

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

        const std::string& quality() const {
            return quality_;
        }

        const std::string& rc() {
            if (rc_.size() != sequence_.size()) create_rc();
            return rc_;
        }

        void trim_sequence(uint32_t begin, uint32_t end);

        friend std::unique_ptr<Read> createRead(uint64_t id, const char* name,
            uint32_t name_length, const char* sequence, uint32_t sequence_length,
            const char* quality, uint32_t quality_length);

        friend bioparser::FastaReader<Read>;
        friend bioparser::FastqReader<Read>;

    private:
        Read(uint64_t id, const char* name, uint32_t name_length, const char* sequence,
            uint32_t sequence_length);
        Read(uint64_t id, const char* name, uint32_t name_length, const char* sequence,
            uint32_t sequence_length, const char* quality, uint32_t quality_length);
        Read(const Read&) = delete;
        const Read& operator=(const Read&) = delete;

        void create_rc();

        uint64_t id_;
        std::string name_;
        std::string sequence_;
        std::string quality_;
        std::string rc_;
};

class ReadInfo;
std::unique_ptr<ReadInfo> createReadInfo(uint64_t id, uint32_t read_length,
    std::vector<uint32_t>& mappings);

std::unique_ptr<ReadInfo> copyReadInfo(std::shared_ptr<ReadInfo> read_info);

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
         * @brief Adds new mappings to coverage_graph_
         */
        void update_coverage_graph(std::vector<uint32_t>& mappings);

        /*!
         * @brief Clears coverage_graph_ and sets begin_, end_ to 0, coverage_graph_.size()
         * respectively
         */
        void reset_coverage_graph();

        /*!
         * @biref Smooths coverage graph with 1D average filter
         */
        void smooth_coverage_graph();

        /*!
         * @brief Correct coverage graph with other graph over overlapping region
         */
        void correct_coverage_graph(uint32_t region_begin, uint32_t region_end,
            std::shared_ptr<ReadInfo> other, uint32_t other_region_begin,
            uint32_t other_region_end, bool rc);

        /*!
         * @brief Locates region in coverage_graph_ with values greater or equal to coverage;
         * updates begin_, end_ and coverage_graph_ accordingly; if there is no such region
         * (with valid coverage and longer than 500), both begin_ and end_ are set to same
         * value and coverage_graph_ is deleted (the read is not valid)
         */
        void find_valid_region(uint32_t coverage);

        /*!
         * @brief Locates pits in coverage_graph_ which ought to indicate that the
         * read is chimeric; if a pit is present, both begin_ and end_ are set to same
         * value and coverage_graph_ is deleted (the read is not valid)
         */
        void find_coverage_pits(double slope_ratio, uint32_t min_slope_width,
            double slope_width_ratio);

        /*!
         * @brief Returns regions of read which ought to be repetitive in the genome
         */
        const std::vector<std::pair<uint32_t, uint32_t>>& coverage_hills() const {
            return coverage_hills_;
        }

        /*!
         * @brief Locates regions in coverage_graph_ which ought to be repetitive in the genome
         * and stores them in coverage_hills_
         */
        void find_coverage_hills(double slope_ratio, uint32_t min_slope_width,
            double slope_width_ratio, double hill_width_ratio, uint32_t dataset_median);

        void find_coverage_hills_simple(uint32_t min_coverage);

        /*!
         * @brief Print coverage_graph_ in csv format to path
         */
        void print_csv(std::string path, uint32_t dataset_median) const;

        friend std::unique_ptr<ReadInfo> createReadInfo(uint64_t id, uint32_t read_length,
            std::vector<uint32_t>& mappings);

        friend std::unique_ptr<ReadInfo> copyReadInfo(std::shared_ptr<ReadInfo> read_info);

    private:
        ReadInfo(uint64_t id, uint32_t read_length, std::vector<uint32_t>& mappings);
        ReadInfo(const ReadInfo&) = default;
        const ReadInfo& operator=(const ReadInfo&) = delete;

        static void coverage_window_add(std::deque<std::pair<int32_t, int32_t>>& window, int32_t value, int32_t position);
        static void coverage_window_update(std::deque<std::pair<int32_t, int32_t>>& window, int32_t position);

        uint64_t id_;
        uint32_t begin_;
        uint32_t end_;
        uint16_t coverage_median_;
        std::vector<uint16_t> coverage_graph_;
        std::vector<std::pair<uint32_t, uint32_t>> coverage_hills_;
};

}
