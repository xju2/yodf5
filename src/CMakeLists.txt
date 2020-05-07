cmake_minimum_required (VERSION 3.10)
project (YODF5 VERSION 0.1)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wfatal-errors)
add_compile_options(-Wno-deprecated)

find_package(HDF5)
set                     (libraries ${libraries} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})


if (NOT HIGHFIVE_DIR)
    message(WARNING "Need HighFive to build program, set to /usr/local")
    set(HIGHFIVE_DIR /usr/local)
endif ()

if (NOT YODA_DIR)
    set(YODA_DIR /usr/local)
endif()

find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# find_package(Boost)

include_directories(
    "/usr/include/mpich-x86_64"
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${HIGHFIVE_DIR}/include
    ${YODA_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)
add_subdirectory(libyodf5)

add_subdirectory(test)
add_subdirectory(utils)


install(EXPORT YODF5Targets
    FILE YODF5Targets.cmake
    DESTINATION lib/cmake/YODF5
)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/YODF5Config.cmake"
    INSTALL_DESTINATION "lib/cmake/example"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/YODF5ConfigVersion.cmake"
  VERSION "${YODF5_VERSION_MAJOR}.${YODF5_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)
# install the configuration file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/YODF5Config.cmake
  DESTINATION lib/cmake/YODF5
)

export(EXPORT YODF5Targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/YODF5Targets.cmake"
)

include(InstallRequiredSystemLibraries)
# set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${YODF5_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${YODF5_VERSION_MINOR}")
include(CPack)