cmake_minimum_required (VERSION 2.8)

# Prevent in-source builds
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(codec2)


#
# Set project version information. This should probably be done via external
# file at some point.
#
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 2)
# Set to patch level is needed, otherwise leave FALSE.
set(PROJECT_VERSION_PATCH FALSE)
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
# Patch level version bumps should not change API/ABI.
set(SOVERSION ${PROJECT_VERSION})
if(PROJECT_VERSION_PATCH)
	set(PROJECT_VERSION "${PROJECT_VERSION}.${PROJECT_VERSION_PATCH}")
endif()

message(STATUS "codec2 version: ${PROJECT_VERSION}")

# Set default C++ flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -O3 -g")

# Build shared libraries by default
set(BUILD_SHARED_LIBS TRUE)


#
# Find some standard headers and functions.
# Should we error out if one is not available?
# Do the results need to be written to a config.h?
#
find_library(LIBM m )
if(LIBM)
	set(CMAKE_REQUIRED_INCLUDES math.h)
	set(CMAKE_REQUIRED_LIBRARIES m)
endif(LIBM)

include(CheckIncludeFiles)
check_include_files("stdlib.h" HAVE_STDLIB_H)
check_include_files("string.h" HAVE_STRING_H)

include(CheckFunctionExists)
check_function_exists(floor HAVE_FLOOR)
check_function_exists(ceil  HAVE_CEIL)
check_function_exists(pow   HAVE_POW)
check_function_exists(sqrt  HAVE_SQRT)
check_function_exists(sin   HAVE_SIN)
check_function_exists(cos   HAVE_COS)
check_function_exists(atan2 HAVE_ATAN2)
check_function_exists(log10 HAVE_LOG10)


#
# codec2 library
#
add_subdirectory(src)

#
# Install example code? Defaults to no.
#
option(INSTALL_EXAMPLES "Install example code." OFF)
if(INSTALL_EXAMPLES)
	install(DIRECTORY octave raw script voicing wav
		DESTINATION share/codec2)
endif()

#
# Cpack configuration to go here...
#
