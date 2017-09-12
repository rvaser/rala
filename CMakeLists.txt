cmake_minimum_required(VERSION 3.2)
project(rala)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

option(rala_build_tests "Build rala google tests" OFF)

add_executable(rala
    src/graph.cpp
    src/main.cpp
    src/overlap.cpp
    src/preprocess.cpp
    src/read.cpp
    src/timer.cpp)

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

target_link_libraries(rala bioparser thread_pool pthread)

if (rala_build_tests)

endif(rala_build_tests)

install(TARGETS rala DESTINATION bin)
