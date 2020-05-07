cmake_minimum_required (VERSION 3.1)
project (YODF5)

find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)

# The version number
set (YODF5_VERSION_MAJOR 0)
set (YODF5_VERSION_MINOR 1)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# source code
add_subdirectory(src)