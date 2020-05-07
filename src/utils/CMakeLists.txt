
add_executable(readH5YODA readH5YODA.cxx)
target_link_libraries(readH5YODA ${libraries} ${YODA_DIR}/lib/libYODA${CMAKE_SHARED_LIBRARY_SUFFIX})


add_executable(readAndWrite readAndWrite.cc)
target_link_libraries(readAndWrite ${libraries} ${YODA_DIR}/lib/libYODA${CMAKE_SHARED_LIBRARY_SUFFIX})


add_executable(pdf_sys pdf_sys.cxx)
target_link_libraries(pdf_sys ${libraries} "/usr/local/lib/libLHAPDF.so")