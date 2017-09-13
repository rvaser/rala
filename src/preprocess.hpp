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
class ReadInfo;
class Overlap;

/*!
 * @brief Removes self and duplicate overlaps, removes overlaps in chimeric/repeat
 * areas, removes contained reads and their overlaps, classifies overlaps;
 * finds contiguos regions in reads which have coverage larger than a predefined
 * threshold (Li 2016)
 */
void preprocessData(std::vector<std::shared_ptr<Read>>& reads,
    std::vector<std::shared_ptr<ReadInfo>>& read_infos,
    std::vector<std::shared_ptr<Overlap>>& overlaps,
    double& dataset_coverage_median, const std::string& reads_path,
    const std::string& overlaps_path, uint32_t overlap_type,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool);

}
