/*!
 * @file preprocess.hpp
 *
 * @brief Preprocessing header file
 */

#include <stdint.h>
#include <string>
#include <memory>
#include <vector>

namespace thread_pool {
    class ThreadPool;
}

namespace rala {

class Read;
class Overlap;

constexpr uint32_t kChunkSize = 1024 * 1024 * 1024; // ~ 1GB

/*
 * @brief Removes self and duplicate overlaps, removes overlaps in chimeric/repeat
 *     areas, removes contained reads and their overlaps, classifies overlaps.
 *     If prefilter is true, finds contiguos regions in reads which have coverage
 *     larger than a predefined threshold (Li 2016)
 */
void preprocessData(std::vector<std::shared_ptr<Read>>& reads, std::vector<std::shared_ptr<Overlap>>& overlaps,
    const std::string& reads_path, const std::string& overlaps_path, uint32_t overlap_type,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool, bool prefilter = true);

/*
 * @brief Chimeric read detection with reference
 */
void findChimericReads(const std::string& reads_path, const std::string& overlaps_path,
    uint32_t overlap_type);

/*
 * @brief Finds overlaps which aren't fully mapped to reference
 */
void findUncontainedReads(const std::string& reads_path, const std::string& overlaps_path,
    uint32_t overlap_type);

/*
 * @brief Creates FASTQ file from FASTA file with dummy quality
 */
void fastaToFastq(const std::string& reads_path);

/*
 * @brief Checks if two numbers are comparable
 */
bool comparable(double a, double b, double eps);

}
