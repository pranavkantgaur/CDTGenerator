add_executable (cdtGenerator rply.cpp cdtGen.cpp)

target_link_libraries (cdtGenerator gmp mpfr CGAL)
