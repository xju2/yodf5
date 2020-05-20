
add_executable(readH5YODA readH5YODA.cxx)
target_link_libraries(readH5YODA ${libraries} ${YODA_DIR}/lib/libYODA${CMAKE_SHARED_LIBRARY_SUFFIX})


add_executable(readAndWrite readAndWrite.cc)
target_link_libraries(readAndWrite ${libraries} ${YODA_DIR}/lib/libYODA${CMAKE_SHARED_LIBRARY_SUFFIX} YODF5)


add_executable(theory_sys theory_sys.cxx)
target_link_libraries(theory_sys ${libraries} YODF5)

install(TARGETS readH5YODA readAndWrite theory_sys 
        DESTINATION bin)