enable_testing ()
include_directories (${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
# add test cpp file
add_executable (CDTGeneratorTest ../src/rply.cpp ../src/cdtGen.cpp testSegmentRecovery.cpp main.cpp)
# link test executable against gtest and gtest-main
target_link_libraries (CDTGeneratorTest gmp mpfr CGAL gtest gtest_main)
add_test (CDTGeneratorTest CDTGeneratorTest)

