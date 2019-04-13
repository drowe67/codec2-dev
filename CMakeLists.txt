#
# Codec2 - Next-Generation Digital Voice for Two-Way Radio
#
# CMake configuration contributed by Richard Shaw (KF5OIM)
# Please report questions, comments, problems, or patches to the freetel
# mailing list: https://lists.sourceforge.net/lists/listinfo/freetel-codec2
#
project(codec2 C)

cmake_minimum_required(VERSION 3.0)

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
set(CODEC2_VERSION_MINOR 8)
# Set to patch level if needed, otherwise leave FALSE.
# Must be positive (non-zero) if set, since 0 == FALSE in CMake.
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

# Set default C flags.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-strict-overflow")

# Check for what C standard is supported.
include(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-std=gnu11" COMPILER_SUPPORTS_GNU11)
CHECK_C_COMPILER_FLAG("-std=gnu99" COMPILER_SUPPORTS_GNU99)

if(COMPILER_SUPPORTS_GNU11)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11")
elseif(COMPILER_SUPPORTS_GNU99)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
else()
   message(SEND_ERROR "Compiler doesn't seem to support at least gnu99, might cause problems" )
endif()

# -fPIC is implied on MinGW...
if(NOT WIN32)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif()

set(CMAKE_C_FLAGS_DEBUG "-g -O2")
set(CMAKE_C_FLAGS_RELEASE "-O3")

#
# Setup Windows/MinGW specifics here.
#
if(MINGW)
    message(STATUS "System is MinGW.")
endif(MINGW)


#
# Find the svn revision if this is a working copy.
#
if(EXISTS ${CMAKE_SOURCE_DIR}/.git)
    find_package(Git)
    if(Git_FOUND)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" describe --always HEAD
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            RESULT_VARIABLE res
            OUTPUT_VARIABLE CODEC2_HASH
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Codec2 current git hash: ${CODEC2_HASH}")
    else()
        message(WARNING "Git not found. Can not determine current commit hash.")
    endif()
endif()


#
# Default options
#
option(BUILD_SHARED_LIBS
    "Build shared library. Set to OFF for static library." ON)
# Unittest should be on for dev builds and off for releases.
if(CMAKE_BUILD_TYPE MATCHES "Release")
    option(UNITTEST "Build unittest binaries." OFF)
else()
    option(UNITTEST "Build unittest binaries." ON)
endif()
option(INSTALL_EXAMPLES "Install example code." OFF)
if(INSTALL_EXAMPLES)
    install(DIRECTORY octave raw script wav
        USE_SOURCE_PERMISSIONS
        DESTINATION ${CMAKE_INSTALL_DATADIR}/codec2)
endif()

# LPCNet needs to be bootstraped because codec2 and freedvlpcnet are 
# cross dependent.
option(LPCNET "Build codec2 with LPCNet support." OFF)
# Setting LPCNET_BUILD_DIR implies LPCNET=ON
if(LPCNET_BUILD_DIR)
    set(LPCNET ON)
endif()


# Math library is automatic on windows
if(POLICY CMP0075)
    cmake_policy(SET CMP0075 NEW)
endif()
if(UNIX)
    set(CMAKE_REQUIRED_LIBRARIES m)
endif(UNIX)

include(CheckIncludeFiles)
check_include_files("stdlib.h" HAVE_STDLIB_H)
check_include_files("string.h" HAVE_STRING_H)

include(CheckSymbolExists)
# Check if _GNU_SOURCE is available.
if (NOT DEFINED _GNU_SOURCE)
  check_symbol_exists(__GNU_LIBRARY__ "features.h" _GNU_SOURCE)

  if (NOT _GNU_SOURCE)
    unset(_GNU_SOURCE CACHE)
    check_symbol_exists(_GNU_SOURCE "features.h" _GNU_SOURCE)
  endif()
endif()

if (_GNU_SOURCE)
    add_definitions(-D_GNU_SOURCE=1)
endif()

check_symbol_exists(floor  math.h HAVE_FLOOR)
check_symbol_exists(ceil   math.h HAVE_CEIL)
check_symbol_exists(pow    math.h HAVE_POW)
check_symbol_exists(sqrt   math.h HAVE_SQRT)
check_symbol_exists(sin    math.h HAVE_SIN)
check_symbol_exists(cos    math.h HAVE_COS)
check_symbol_exists(atan2  math.h HAVE_ATAN2)
check_symbol_exists(log10  math.h HAVE_LOG10)
check_symbol_exists(round  math.h HAVE_ROUND)
check_symbol_exists(getopt math.h HAVE_GETOPT)

configure_file ("${PROJECT_SOURCE_DIR}/cmake/config.h.in"
                "${PROJECT_BINARY_DIR}/config.h" )
# Output path is such that #include <codec2/version.h> in codec2.h works
set(CODEC2_VERSION_PATH "${PROJECT_BINARY_DIR}/codec2")
configure_file ("${PROJECT_SOURCE_DIR}/cmake/version.h.in"
                "${CODEC2_VERSION_PATH}/version.h" )
include_directories(${PROJECT_BINARY_DIR})

# CMake Package setup
#include(CMakePackageConfigHelpers)
#configure_package_config_file(cmake/codec2-config.cmake.in
#    ${CMAKE_CURRENT_BINARY_DIR}/codec2-config.cmake
#    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/codec2
#    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
#)

#
# Find lpcnet library
#
if(LPCNET)
    if(LPCNET_BUILD_DIR)
        message(STATUS "Using LPCNet from build directory...")
        message(STATUS "  LPCNet Build dir: ${LPCNET_BUILD_DIR}")
        load_cache(${LPCNET_BUILD_DIR} LPCNet_SOURCE_DIR)
        find_path(LPCNET_INCLUDE_DIR lpcnet_freedv.h
            PATHS ${LPCNet_SOURCE_DIR}
            PATH_SUFFIXES src
            NO_DEFAULT_PATH)
        message(STATUS "  Source include dir: ${LPCNET_INCLUDE_DIR}")
        find_library(LPCNET_LIBRARY NAME lpcnetfreedv
            PATHS ${LPCNET_BUILD_DIR}
            PATH_SUFFIXES src
            NO_DEFAULT_PATH)
        message(STATUS "  Library: ${LPCNET_LIBRARY}")
        if(LPCNET_INCLUDE_DIR AND LPCNET_LIBRARY)
            add_definitions("-D__LPCNET__")
            include_directories(${LPCNET_INCLUDE_DIR})
            set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES} ${LPCNET_LIBRARY}")
        else()
            message(FATAL_ERROR "LPCNet include/library not found.")
        endif()
    else()
        find_package(lpcnet REQUIRED)
        if(lpcnet_FOUND)
            add_definitions("-D__LPCNET__")
        else()
            message(FATAL_ERROR "lpcnetfreedv library not found.")
        endif()
    endif()
endif()

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
    find_path(SPEEXDSP_INCLUDE_DIR speex/speex_preprocess.h)
    find_library(SPEEXDSP_LIBRARY speexdsp)
    message(STATUS "  Speex DSP headers: ${SPEEXDSP_INCLUDE_DIR}")
    message(STATUS "  Speex DSP library: ${SPEEXDSP_LIBRARY}")
    if(NOT SPEEXDSP_INCLUDE_DIR AND NOT SPEEXDSP_LIBRARY)
        message(FATAL_ERROR "Speex DSP library not found!")
    endif()

    #
    # Samplerate Library
    #
    message(STATUS "Looking for samplerate...")
    find_library(LIBSAMPLERATE samplerate)
    find_path(LIBSAMPLERATE_INCLUDE_DIR samplerate.h)
    message(STATUS "  samplerate headers: ${LIBSAMPLERATE_INCLUDE_DIR}")
    message(STATUS "  samplerate library: ${LIBSAMPLERATE}")
    if(LIBSAMPLERATE AND LIBSAMPLERATE_INCLUDE_DIR)
        set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES} ${LIBSAMPLERATE}")
        include_directories(${LIBSAMPLERATE_INCLUDE_DIR})
    else(LIBSTAMPLERATE AND LIBSAMPLERATE_INCLUDE_DIR)
        message(FATAL_ERROR "samplerate library not found.
On Linux systems try installing:
    samplerate-devel  (RPM based systems)
    libsamplerate-dev (DEB based systems)")
    endif(LIBSAMPLERATE AND LIBSAMPLERATE_INCLUDE_DIR)

    add_subdirectory(unittest)
    add_subdirectory(misc)
endif(UNITTEST)

message(STATUS "Build type is: " ${CMAKE_BUILD_TYPE})
string(TOUPPER ${CMAKE_BUILD_TYPE} _FLAGS)
if(_FLAGS STREQUAL "NONE")
    message(STATUS "Compiler Flags: " ${CMAKE_C_FLAGS})
else()
    message(STATUS "Compiler Flags: " ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${_FLAGS}})
endif()
message(STATUS "Libraries linked: " ${CMAKE_REQUIRED_LIBRARIES})

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

##################################################################
# Tests
##################################################################

if(UNITTEST)
    include(CTest)
    enable_testing()

    add_test(NAME test_FDMDV_modem_octave_port
             COMMAND sh -c "$<TARGET_FILE:tfdmdv> && octave --no-gui -qf ${CMAKE_CURRENT_SOURCE_DIR}/octave/tfdmdv.m"
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/octave)
             set_tests_properties(test_FDMDV_modem_octave_port PROPERTIES
             PASS_REGULAR_EXPRESSION "fails: 0")

    add_test(NAME test_OFDM_modem_octave_port
             COMMAND sh -c "PATH_TO_TOFDM=${CMAKE_CURRENT_BINARY_DIR}/unittest/tofdm octave --no-gui -qf ${CMAKE_CURRENT_SOURCE_DIR}/octave/tofdm.m"
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/octave)
             set_tests_properties(test_OFDM_modem_octave_port PROPERTIES
             PASS_REGULAR_EXPRESSION "fails: 0")

    add_test(NAME test_OFDM_modem_octave_port_Nc_31
             COMMAND sh -c "NC=31 PATH_TO_TOFDM=${CMAKE_CURRENT_BINARY_DIR}/unittest/tofdm octave --no-gui -qf ${CMAKE_CURRENT_SOURCE_DIR}/octave/tofdm.m"
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/octave)
             set_tests_properties(test_OFDM_modem_octave_port PROPERTIES
             PASS_REGULAR_EXPRESSION "fails: 0")

   add_test(NAME test_OFDM_modem_AWGN_BER
             COMMAND sh -c "$<TARGET_FILE:ofdm_mod> --in /dev/zero --ldpc 1 --testframes 60 --txbpf | $<TARGET_FILE:cohpsk_ch> - - -20 --Fs 8000 -f -10 | $<TARGET_FILE:ofdm_demod> --out /dev/null --testframes --ldpc 1 --verbose 1"
             )

    add_test(NAME test_horus_binary
             COMMAND sh -c "sox ${CMAKE_CURRENT_SOURCE_DIR}/wav/horus_binary_ebno_4.5db.wav -r 48000 -t raw - | $<TARGET_FILE:horus_demod> -m binary - -"
             )
             set_tests_properties(test_horus_binary PROPERTIES
             PASS_REGULAR_EXPRESSION "1C9A9545")

    # tests ldpc_enc/ldpc_noise/ldpc_dec
    add_test(NAME test_ldpc_enc_dec
             COMMAND sh -c "$<TARGET_FILE:ldpc_enc> /dev/zero - --sd --code HRA_112_112 --testframes 200 | $<TARGET_FILE:ldpc_noise> - - 0.5 | $<TARGET_FILE:ldpc_dec> - /dev/null --code HRA_112_112 --sd --testframes"
             )

    add_test(NAME test_ldpc_enc_dec_HRAb_396_504
             COMMAND sh -c "$<TARGET_FILE:ldpc_enc> /dev/zero - --sd --code HRAb_396_504 --testframes 200 | $<TARGET_FILE:ldpc_noise> - - -2.0 | $<TARGET_FILE:ldpc_dec> - /dev/null --code HRAb_396_504 --sd --testframes"
             )

    # some more LDPC tests, placed in a script for convenience
    add_test(NAME test_ldpc_check
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/unittest/ldpc_check
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unittest
             )
             set_tests_properties(test_ldpc_check PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

    # basic test of C utils
    add_test(NAME test_FSK_modem_mod_demod
             COMMAND sh -c "$<TARGET_FILE:fsk_get_test_bits> - 10000 | $<TARGET_FILE:fsk_mod> 2 8000 100 1200 100 - - | $<TARGET_FILE:fsk_demod> -l 2 8000 100 - - | $<TARGET_FILE:fsk_put_test_bits> -"
             )
             set_tests_properties(test_FSK_modem_mod_demod  PROPERTIES PASS_REGULAR_EXPRESSION "bits tested 9800, bit errors 0")

    add_test(NAME test_FreeDV_API_700D_AWGN_BER
             COMMAND sh -c "dd bs=1280 count=120 if=/dev/zero | $<TARGET_FILE:freedv_tx> 700D - - --testframes | $<TARGET_FILE:cohpsk_ch> - - -20 --Fs 8000 -f -10 | $<TARGET_FILE:freedv_rx> 700D - /dev/null --testframes --discard"
             )

    add_test(NAME test_FreeDV_API_700D_AWGN_BER_USECOMPLEX
             COMMAND sh -c "dd bs=1280 count=120 if=/dev/zero | $<TARGET_FILE:freedv_tx> 700D - - --testframes | $<TARGET_FILE:cohpsk_ch> - - -20 --Fs 8000 -f -10 | $<TARGET_FILE:freedv_rx> 700D - /dev/null --testframes --discard --usecomplex"
             )
           
    add_test(NAME test_FreeDV_API_700D_AWGN_BER_Interleave
             COMMAND sh -c "dd bs=1280 count=120 if=/dev/zero | $<TARGET_FILE:freedv_tx> 700D - - --interleave 8 --testframes | $<TARGET_FILE:cohpsk_ch> - - -20 --Fs 8000 -f -10 | $<TARGET_FILE:freedv_rx> 700D - /dev/null --interleave 8 --testframes"
             )

endif(UNITTEST)
