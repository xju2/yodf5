include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )

file(GLOB SrcFiles *.c*)
add_library(YODF5 SHARED ${SrcFiles})
target_include_directories(YODF5 
        INTERFACE 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        )
target_link_libraries(YODF5 LHAPDF)

install(TARGETS YODF5 
        DESTINATION lib
        EXPORT YODF5Targets)

file(GLOB HeaderFiles include/yodf5/*.hpp)
install(FILES ${HeaderFiles}
    DESTINATION include)