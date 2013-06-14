#
# Codec2 - Next-Generation Digital Voice for Two-Way Radio
#
# CMake configuration contributed by Richard Shaw (KF5OIM)
# Please report questions, comments, problems, or patches to the freetel
# mailing list: https://lists.sourceforge.net/lists/listinfo/freetel-codec2
#

cmake_minimum_required(VERSION 2.8)

#
# Prevent in-source builds
# If an in-source build is attempted, you will still need to clean up a few
# files manually.
#
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(codec2)


#
# Set project version information. This should probably be done via external
# file at some point.
#
set(CODEC2_VERSION_MAJOR 0)
set(CODEC2_VERSION_MINOR 2)
# Set to patch level is needed, otherwise leave FALSE.
set(CODEC2_VERSION_PATCH FALSE)
set(CODEC2_VERSION "${CODEC2_VERSION_MAJOR}.${CODEC2_VERSION_MINOR}")
# Patch level version bumps should not change API/ABI.
set(SOVERSION "${CODEC2_VERSION_MAJOR}.${CODEC2_VERSION_MINOR}")
if(CODEC2_VERSION_PATCH)
    set(PROJECT_VERSION "${CODEC2_VERSION}.${CODEC2_VERSION_PATCH}")
endif(CODEC2_VERSION_PATCH)
message(STATUS "codec2 version: ${CODEC2_VERSION}")

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")

#
# Setup Windows/MinGW specifics here.
#
if(MINGW)
    message(STATUS "System is MinGW.")
    foreach(RUNTIME
        libgcc_s_sjlj-1.dll
        libstdc++-6.dll)
        message(STATUS "Checking for ${RUNTIME}")
        find_library(${RUNTIME}_LIB ${RUNTIME})
        message(STATUS "runtime found: ${${RUNTIME}_LIB}")
        list(APPEND CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${${RUNTIME}_LIB})
    endforeach()
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

# Set default C++ flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall")

#
# Default options
#a
option(BUILD_SHARED_LIBS
    "Build shared library. Set to OFF for static library." ON)
option(UNITTEST "Build unittest binaries." OFF)
option(INSTALL_EXAMPLES "Install example code." OFF)
if(INSTALL_EXAMPLES)
    install(DIRECTORY octave raw script voicing wav
        DESTINATION share/codec2)
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

#
# codec2 library
#
add_subdirectory(src)

if(UNITTEST)
    # Pthread Library
    find_package(Threads REQUIRED)
    message(STATUS "Threads library flags: ${CMAKE_THREAD_LIBS_INIT}")
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
    include(InstallRequiredSystemLibraries)
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Next-Generation Digital Voice for Two-Way Radio")
    set(CPACK_PACKAGE_VENDOR "CMake")
    set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
    set(CPACK_PACKAGE_VERSION_MAJOR ${CODEC2_VERSION_MAJOR})
    set(CPACK_PACKAGE_VERSION_MINOR ${CODEC2_VERSION_MINOR})
    if(CODEC2_VERSION_PATCH)
        set(CPACK_PACKAGE_VERSION_PATCH ${CODEC2_VERSION_PATCH})
    endif(CODEC2_VERSION_PATCH)
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "Codec2")
    # There is a bug in NSI that does not handle full unix paths properly. Make
    # sure there is at least one set of four (4) backlasshes.
    #set(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
    #set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
    set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
    set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\codec2.org")
    set(CPACK_NSIS_MODIFY_PATH ON)
    include(CPack)
endif(WIN32)
