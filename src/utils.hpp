/*!
 * @file utils.hpp
 *
 * @brief Utils header file
 */

#include <stdint.h>
#include <string>

namespace rala {

/*
 * @brief Checks if two numbers are similar
 */
bool isSimilar(double a, double b, double eps);

/*
 * @brief Chimeric read detection with reference
 */
void findChimericReads(const std::string& reads_path, const std::string& overlaps_path, uint32_t overlap_type);

/*
 * @brief Finds overlaps which aren't fully mapped to reference
 */
void findUnusedOverlaps(const std::string& reads_path, const std::string& overlaps_path, uint32_t overlap_type);

}
