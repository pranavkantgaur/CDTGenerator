add_executable (cdtGenerator rply.cpp cdtGen.cpp main.cpp)

IF (UNIX)
	target_link_libraries (cdtGenerator gmp mpfr boost_thread CGAL)
ELSEIF (MSVC)
	IF (${MSVC_VERSION} MATCHES 1700) # visual studio 2012
		target_link_libraries (cdtGenerator gmp mpfr boost_thread-vc110-mt-gd-1_58)
	ENDIF()
ENDIF()	
	
