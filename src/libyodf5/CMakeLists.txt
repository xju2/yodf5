include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )

file(GLOB SrcFiles *.c*)
add_library(YODF5 SHARED ${SrcFiles})