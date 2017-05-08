/*!
 * @file timer.cpp
 *
 * @brief Timer class source file
 */

#include <stdio.h>

#include "timer.hpp"

namespace rala {

Timer::Timer()
    : paused_(false), time_(0), timeval_() {
}

void Timer::start() {
    gettimeofday(&timeval_, nullptr);
    paused_ = false;
}

void Timer::stop() {
    if (paused_) {
        return;
    }

    timeval stop;
    gettimeofday(&stop, nullptr);
    time_ += ((stop.tv_sec - timeval_.tv_sec) * 1000000L + stop.tv_usec)
        - timeval_.tv_usec;
    paused_ = true;
}

void Timer::reset() {
    gettimeofday(&timeval_, nullptr);
    time_ = 0;
    paused_ = false;
}

void Timer::print(const char* message) const {
    fprintf(stderr, "%s %.5lf s\n", message, time_ / (double) 1000000);
}

}
