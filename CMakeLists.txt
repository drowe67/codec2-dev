#
# Codec2 - Next-Generation Digital Voice for Two-Way Radio
#
# CMake configuration contributed by Richard Shaw (KF5OIM)
# Please report questions, comments, problems, or patches to the freetel
# mailing list: https://lists.sourceforge.net/lists/listinfo/freetel-codec2
#
project(codec2 C)

cmake_minimum_required(VERSION 2.8)

include(GNUInstallDirs)
mark_as_advanced(CLEAR
    CMAKE_INSTALL_BINDIR
    CMAKE_INSTALL_INCLUDEDIR
    CMAKE_INSTALL_LIBDIR
)

#
# Prevent in-source builds
# If an in-source build is attempted, you will still need to clean up a few
# files manually.
#
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds in ${CMAKE_BINARY_DIR} are not "
   "allowed, please remove ./CMakeCache.txt and ./CMakeFiles/, create a "
   "separate build directory and run cmake from there.")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")


#
# Set project version information. This should probably be done via external
# file at some point.
#
set(CODEC2_VERSION_MAJOR 0)
set(CODEC2_VERSION_MINOR 6)
# Set to patch level if needed, otherwise leave FALSE.
set(CODEC2_VERSION_PATCH FALSE)
set(CODEC2_VERSION "${CODEC2_VERSION_MAJOR}.${CODEC2_VERSION_MINOR}")
# Patch level version bumps should not change API/ABI.
set(SOVERSION "${CODEC2_VERSION_MAJOR}.${CODEC2_VERSION_MINOR}")
if(CODEC2_VERSION_PATCH)
    set(CODEC2_VERSION "${CODEC2_VERSION}.${CODEC2_VERSION_PATCH}")
endif()
message(STATUS "codec2 version: ${CODEC2_VERSION}")

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Set default C++ flags.
include(CheckCCompilerFlag)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -g -O3")
CHECK_C_COMPILER_FLAG("-std=gnu11" COMPILER_SUPPORTS_GNU11)
if(COMPILER_SUPPORTS_GNU11)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
endif()
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS}")

# -fPIC is implied on MinGW...
if(NOT WIN32)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif()

message(STATUS "Build type is: " ${CMAKE_BUILD_TYPE})
string(TOUPPER ${CMAKE_BUILD_TYPE} _FLAGS)
if(NOT _FLAGS STREQUAL "NONE")
    set(BUILD_FLAGS "${CMAKE_C_FLAGS_${_FLAGS}}")
endif()
message(STATUS "Compiler Flags: " ${BUILD_FLAGS})

#
# Setup Windows/MinGW specifics here.
#
if(MINGW)
    message(STATUS "System is MinGW.")
endif(MINGW)


#
# Find the svn revision if this is a working copy.
# WORK IN PROGRESS
# Works ok if it is a working copy but errors out if not.
#
#find_package(Subversion)
#if(Subversion_FOUND)
#   Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} CODEC2)
#   message(STATUS "codec2 svn revision: ${CODEC2_WC_REVISION}")
#else(SUBVERSION_FOUND)
#   message(WARNING "Subversion not found. Can not determine svn revision.")
#endif(SUBVERSION_FOUND)


#
# Default options
#
option(BUILD_SHARED_LIBS
    "Build shared library. Set to OFF for static library." ON)
# Unittest should be on for dev builds and off for releases.
#if(CMAKE_BUILD_TYPE MATCHES "Release")
    option(UNITTEST "Build unittest binaries." OFF)
#else()
    #option(UNITTEST "Build unittest binaries." ON)
#endif()
option(INSTALL_EXAMPLES "Install example code." OFF)
if(INSTALL_EXAMPLES)
    install(DIRECTORY octave raw script wav
        USE_SOURCE_PERMISSIONS
        DESTINATION ${CMAKE_INSTALL_DATADIR}/codec2)
endif()


# Math library is automatic on windows
if(UNIX)
    set(CMAKE_REQUIRED_INCLUDES math.h)
    set(CMAKE_REQUIRED_LIBRARIES m)
endif(UNIX)

include(CheckIncludeFiles)
check_include_files("stdlib.h" HAVE_STDLIB_H)
check_include_files("string.h" HAVE_STRING_H)

include(CheckFunctionExists)
check_function_exists(floor  HAVE_FLOOR)
check_function_exists(ceil   HAVE_CEIL)
check_function_exists(pow    HAVE_POW)
check_function_exists(sqrt   HAVE_SQRT)
check_function_exists(sin    HAVE_SIN)
check_function_exists(cos    HAVE_COS)
check_function_exists(atan2  HAVE_ATAN2)
check_function_exists(log10  HAVE_LOG10)
check_function_exists(round  HAVE_ROUND)
check_function_exists(getopt HAVE_GETOPT)

configure_file ("${PROJECT_SOURCE_DIR}/cmake/config.h.in"
                "${PROJECT_BINARY_DIR}/config.h" )
include_directories(${PROJECT_BINARY_DIR})

# CMake Package setup
#include(CMakePackageConfigHelpers)
#configure_package_config_file(cmake/codec2-config.cmake.in
#    ${CMAKE_CURRENT_BINARY_DIR}/codec2-config.cmake
#    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/codec2
#    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
#)

#
# codec2 library
#
add_subdirectory(src)

if(UNITTEST)
    # Pthread Library
    find_package(Threads REQUIRED)
    message(STATUS "Threads library flags: ${CMAKE_THREAD_LIBS_INIT}")

    #
    # Find speex library
    #
    message(STATUS "Looking for Speex DSP library.")
    find_path(SPEEXDSP_INCLUDE_DIR speex/speex.h)
    find_library(SPEEXDSP_LIBRARY speexdsp)
    message(STATUS "  Speex DSP headers: ${SPEEXDSP_INCLUDE_DIR}")
    message(STATUS "  Speex DSP library: ${SPEEXDSP_LIBRARY}")
    if(NOT SPEEXDSP_INCLUDE_DIR AND NOT SPEEXDSP_LIBRARY)
        message(FATAL_ERROR "Speex DSP library not found!")
    endif()

    add_subdirectory(unittest)
endif(UNITTEST)

#
# Cpack NSIS installer configuration for Windows.
# See: http://nsis.sourceforge.net/Download
#
# *nix systems should use "make install" and/or appropriate
# distribution packaging tools.
#
if(WIN32)
    # Detect if we're doing a 32-bit or 64-bit windows build.
    if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
        set(CMAKE_CL_64 TRUE)
    endif()
    configure_file(cmake/GetDependencies.cmake.in cmake/GetDependencies.cmake
        @ONLY
    )
    install(SCRIPT ${CMAKE_BINARY_DIR}/cmake/GetDependencies.cmake)
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Next-Generation Digital Voice for Two-Way Radio")
    set(CPACK_PACKAGE_VENDOR "CMake")
    set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
    set(CPACK_PACKAGE_VERSION_MAJOR ${CODEC2_VERSION_MAJOR})
    set(CPACK_PACKAGE_VERSION_MINOR ${CODEC2_VERSION_MINOR})
    if(CODEC2_VERSION_PATCH)
        set(CPACK_PACKAGE_VERSION_PATCH ${CODEC2_VERSION_PATCH})
    else()
        set(CPACK_PACKAGE_VERSION_PATCH 0)
    endif()
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "Codec2")
    set(CPACK_CREATE_DESKTOP_LINKS "")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
    set(CPACK_NSIS_URL_INFO_ABOUT "http://codec2.org")
    set(CPACK_NSIS_MODIFY_PATH ON)
    include(CPack)
endif(WIN32)

########################################################################
# Create Pkg Config File
########################################################################
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/codec2.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/codec2.pc
    @ONLY
)

install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/codec2.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    COMPONENT "codec2_devel"
)
