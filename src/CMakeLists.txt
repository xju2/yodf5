set(CMAKE_CXX_STANDARD 14)

add_compile_options(-Wfatal-errors)
add_compile_options(-Wno-deprecated)

find_package(HDF5)
set                     (libraries ${libraries} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})


if (NOT HIGHFIVE_DIR)
    message(FATAL_ERROR "Need HighFive to build program --- check -DHIGHFIVE_DIR")
endif ()

find_package(xtensor)
# find_package(Boost)

include_directories(
    "/usr/include/mpich-x86_64"
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${HIGHFIVE_DIR}/include
    # ${Boost_INCLUDE_DIR}
)
add_subdirectory(test)