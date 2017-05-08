/*!
 * @file timer.cpp
 *
 * @brief Timer class source file
 */

#pragma once

#include <stdint.h>
#include <future>
#include <sys/time.h>

namespace rala {

class Timer {
public:

    Timer();

    /*!
     * @brief Records the current time into timeval_ and unpauses
     * the timer.
     */
    void start();

    /*!
     * @brief Subtracts the current time from the time in timeval_
     * and adds the difference to time_ if the timer is not paused.
     */
    void stop();

    /*!
     * @brief Resets the variable time_ and starts the timer.
     */
    void reset();

    /*!
     * @brief Prints to stderr the elapsed time in seconds in following
     * format: message time (s).
     */
    void print(const char* message) const;

private:

    bool paused_;
    uint64_t time_;
    timeval timeval_;
};

}
