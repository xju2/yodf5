add_executable(test_xtensor test_xtensor.cxx)
target_link_libraries(test_xtensor ${libraries})

add_executable(test_pdf_sys test_pdf_sys.cxx)
target_link_libraries(test_pdf_sys ${libraries} "/usr/local/lib/libLHAPDF.so")