set(CMAKE_CXX_STANDARD 20)

add_compile_options(-Wfatal-errors)
add_compile_options(-Wno-deprecated)

find_package(HDF5)
set                     (libraries ${libraries} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})


if (NOT HIGHFIVE_DIR)
    message(WARNING "Need HighFive to build program --- check -DHIGHFIVE_DIR")
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