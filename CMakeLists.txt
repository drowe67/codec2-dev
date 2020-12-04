#
# Codec2 - Next-Generation Digital Voice for Two-Way Radio
#
# CMake configuration contributed by Richard Shaw (KF5OIM)
# Please report questions, comments, problems, or patches to the freetel
# mailing list: https://lists.sourceforge.net/lists/listinfo/freetel-codec2
#
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "Minimum OS X deployment version")

project(codec2 C)

cmake_minimum_required(VERSION 3.0)

# Set policies here, probably should move to cmake dir.
if(POLICY CMP0075)
    cmake_policy(SET CMP0075 NEW)
endif()
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()

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
set(CODEC2_VERSION_MINOR 9)
# Set to patch level if needed, otherwise leave FALSE.
# Must be positive (non-zero) if set, since 0 == FALSE in CMake.
set(CODEC2_VERSION_PATCH 2)
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

# Build universal ARM64 and x86_64 binaries on Mac.
if(BUILD_OSX_UNIVERSAL)
set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64")
endif(BUILD_OSX_UNIVERSAL)

#
# Find the git hash if this is a working copy.
#
if(EXISTS ${CMAKE_SOURCE_DIR}/.git)
    find_package(Git QUIET)
    if(Git_FOUND)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" describe --always HEAD
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            RESULT_VARIABLE res
            OUTPUT_VARIABLE FREEDV_HASH
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "freedv-gui current git hash: ${FREEDV_HASH}")
        add_definitions(-DGIT_HASH="${FREEDV_HASH}")
    else()
        message(WARNING "Git not found. Can not determine current commit hash.")
        add_definitions(-DGIT_HASH="Unknown")
    endif()
else()
        add_definitions(-DGIT_HASH="None")
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

set(CMAKE_C_FLAGS_DEBUG "-g -O2 -DDUMP")
set(CMAKE_C_FLAGS_RELEASE "-O3")

#
# Setup Windows/MinGW specifics here.
#
if(MINGW)
    message(STATUS "System is MinGW.")
endif(MINGW)


#
# Find the git hash if this is a working copy.
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
set(LPCNET_BUILD_DIR FALSE CACHE PATH "Location of lpcnet build tree.")
# Setting LPCNET_BUILD_DIR implies LPCNET=ON
if(LPCNET_BUILD_DIR)
    set(LPCNET ON)
endif()



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

if(UNIX)
    set(CMAKE_REQUIRED_LIBRARIES m)
endif()

check_symbol_exists(floor  math.h   HAVE_FLOOR)
check_symbol_exists(ceil   math.h   HAVE_CEIL)
check_symbol_exists(pow    math.h   HAVE_POW)
check_symbol_exists(sqrt   math.h   HAVE_SQRT)
check_symbol_exists(sin    math.h   HAVE_SIN)
check_symbol_exists(cos    math.h   HAVE_COS)
check_symbol_exists(atan2  math.h   HAVE_ATAN2)
check_symbol_exists(log10  math.h   HAVE_LOG10)
check_symbol_exists(round  math.h   HAVE_ROUND)
check_symbol_exists(getopt getopt.h HAVE_GETOPT)

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
        find_package(lpcnetfreedv REQUIRED
            PATHS ${LPCNET_BUILD_DIR}
            NO_DEFAULT_PATH
            CONFIGS lpcnetfreedv.cmake
        )
        if(lpcnetfreedv_FOUND)
            message(STATUS "liblpcnetfreedv found in build tree.")
            add_definitions("-D__LPCNET__")
        else()
            message(FATAL_ERROR "LPCNet include/library not found in build tree.")
        endif()
    else()
        find_package(lpcnetfreedv REQUIRED)
        if(lpcnetfreedv_FOUND)
            add_definitions("-D__LPCNET__")
            message(STATUS "liblpcnetfreedv found.")
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
    set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
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
    set(CPACK_NSIS_URL_INFO_ABOUT "http://rowetel.com/codec2.html")
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

    add_test(NAME test_freedv_get_hash
             COMMAND sh -c "${CMAKE_CURRENT_BINARY_DIR}/unittest/thash")

    add_test(NAME test_CML_ldpcut
             COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/octave; SHORT_VERSION_FOR_CTEST=1 octave --no-gui -qf ldpcut.m")
             set_tests_properties(test_CML_ldpcut PROPERTIES PASS_REGULAR_EXPRESSION "Nerr: 0")

    add_test(NAME test_codec2_700c_octave_port
             COMMAND sh -c "
               cd ${CMAKE_CURRENT_BINARY_DIR}/src;
               ./c2sim ${CMAKE_CURRENT_SOURCE_DIR}/raw/cq_ref.raw --phase0 --postfilter --dump cq_ref --lpc 10 --dump_pitch_e cq_ref_pitche.txt;
               cd ${CMAKE_CURRENT_BINARY_DIR}/unittest; ./tnewamp1 ${CMAKE_CURRENT_SOURCE_DIR}/raw/cq_ref.raw;
               cd ${CMAKE_CURRENT_SOURCE_DIR}/octave;
               DISPLAY=\"\" octave-cli -qf --eval 'tnewamp1(\"${CMAKE_CURRENT_BINARY_DIR}/src/cq_ref\", \"${CMAKE_CURRENT_BINARY_DIR}/unittest\")'")
             set_tests_properties(test_codec2_700c_octave_port PROPERTIES PASS_REGULAR_EXPRESSION "fails: 0")

    add_test(NAME test_FDMDV_modem_octave_port
             COMMAND sh -c "$<TARGET_FILE:tfdmdv> && DISPLAY=\"\" octave-cli --no-gui -qf ${CMAKE_CURRENT_SOURCE_DIR}/octave/tfdmdv.m"
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/octave)
             set_tests_properties(test_FDMDV_modem_octave_port PROPERTIES
             PASS_REGULAR_EXPRESSION "fails: 0")

    add_test(NAME test_COHPSK_modem_octave_port
             COMMAND sh -c "$<TARGET_FILE:tcohpsk> && DISPLAY=\"\" octave-cli --no-gui -qf ${CMAKE_CURRENT_SOURCE_DIR}/octave/tcohpsk.m"
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/octave)
             set_tests_properties(test_COHPSK_modem_octave_port PROPERTIES
             PASS_REGULAR_EXPRESSION "fails: 0")

    add_test(NAME test_COHPSK_modem_AWGN_BER
             COMMAND sh -c "$<TARGET_FILE:cohpsk_get_test_bits> - 5600 | $<TARGET_FILE:cohpsk_mod> - - | $<TARGET_FILE:cohpsk_ch> - - -30  | $<TARGET_FILE:cohpsk_demod> - - | $<TARGET_FILE:cohpsk_put_test_bits> -"
             )
           
    add_test(NAME test_COHPSK_modem_freq_offset
             COMMAND sh -c "set -x; $<TARGET_FILE:cohpsk_get_test_bits> - 5600 | $<TARGET_FILE:cohpsk_mod> - - | $<TARGET_FILE:cohpsk_ch> - - -40 -f -20 | $<TARGET_FILE:cohpsk_demod> -v - - 2>log.txt | $<TARGET_FILE:cohpsk_put_test_bits> - ; ! grep 'lost sync' log.txt"
             )
           
    # -------------------------------------------------------------------------
    #                                OFDM Modem
    # -------------------------------------------------------------------------

    add_test(NAME test_OFDM_qam16
             COMMAND sh -c "${CMAKE_CURRENT_BINARY_DIR}/unittest/tqam16")

    add_test(NAME test_OFDM_modem_octave_port
             COMMAND sh -c "PATH_TO_TOFDM=${CMAKE_CURRENT_BINARY_DIR}/unittest/tofdm DISPLAY=\"\" octave-cli --no-gui -qf ${CMAKE_CURRENT_SOURCE_DIR}/octave/tofdm.m"
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/octave)
             set_tests_properties(test_OFDM_modem_octave_port PROPERTIES
             PASS_REGULAR_EXPRESSION "fails: 0")

    add_test(NAME test_OFDM_modem_octave_port_Nc_31
             COMMAND sh -c "NC=31 PATH_TO_TOFDM=${CMAKE_CURRENT_BINARY_DIR}/unittest/tofdm DISPLAY=\"\" octave-cli --no-gui -qf ${CMAKE_CURRENT_SOURCE_DIR}/octave/tofdm.m"
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/octave)
             set_tests_properties(test_OFDM_modem_octave_port_Nc_31 PROPERTIES
             PASS_REGULAR_EXPRESSION "fails: 0")

    # noise free uncoded 700D test, including reading and writing payload bits
    add_test(NAME test_OFDM_modem_700D
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ofdm_get_test_bits - | 
                            ./ofdm_mod | 
                            ./ofdm_demod  --testframes > /dev/null")
           
    # noise free coded 700D test, including reading and writing payload bits
    add_test(NAME test_OFDM_modem_700D_ldpc
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ofdm_get_test_bits - --ldpc | 
                            ./ofdm_mod --ldpc | 
                            ./ofdm_demod --ldpc --testframes > /dev/null")
           
    # noise free 2020 test, including reading and writing payload bits.  fsk_*_test_bits used as it does it's own frame sync
    add_test(NAME test_OFDM_modem_2020_ldpc
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./fsk_get_test_bits - 5000 | 
                            ./ofdm_mod --ldpc --mode 2020 -p 312 | 
                            ./ofdm_demod  --ldpc --mode 2020 -p 312 | 
                            ./fsk_put_test_bits - -q")

    add_test(NAME test_OFDM_modem_AWGN_BER
             COMMAND sh -c "$<TARGET_FILE:ofdm_mod> --in /dev/zero --ldpc --testframes 60 --txbpf | $<TARGET_FILE:cohpsk_ch> - - -20 --Fs 8000 -f -50 | $<TARGET_FILE:ofdm_demod> --out /dev/null --testframes --ldpc --verbose 1"
             )

    configure_file(unittest/ofdm_fade.sh.in unittest/ofdm_fade.sh)
    add_test(NAME test_OFDM_modem_fading_BER
             COMMAND ${CMAKE_CURRENT_BINARY_DIR}/unittest/ofdm_fade.sh
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/unittest
             )
           
    add_test(NAME test_OFDM_modem_phase_est_bw
             COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/unittest;
                            PATH=$PATH:${CMAKE_CURRENT_BINARY_DIR}/src ./ofdm_phase_est_bw.sh")
           
    add_test(NAME test_OFDM_modem_fading_DPSK_BER
             COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/unittest;
                            PATH=$PATH:${CMAKE_CURRENT_BINARY_DIR}/src ./ofdm_fade_dpsk.sh")

    add_test(NAME test_OFDM_modem_time_sync_700D
             COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/unittest;
                            PATH=$PATH:${CMAKE_CURRENT_BINARY_DIR}/src ./ofdm_time_sync.sh 700D")
           
if(LPCNET)
    add_test(NAME test_OFDM_modem_time_sync_2020
             COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/unittest;
                            PATH=$PATH:${CMAKE_CURRENT_BINARY_DIR}/src ./ofdm_time_sync.sh 2020")
endif()

    # -------------------------------------------------------------------------
    #                                OFDM Data modes
    # -------------------------------------------------------------------------

    # C Tx, Octave Rx
    add_test(NAME test_OFDM_modem_datac1_octave
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR};
                            ./src/ofdm_mod --mode datac1 --in /dev/zero --testframes 5 --verbose 1 > test.raw;
                            cd ${CMAKE_CURRENT_SOURCE_DIR}/octave;
                            DISPLAY=\"\" octave-cli -qf --eval 'ofdm_rx(\"${CMAKE_CURRENT_BINARY_DIR}/test.raw\",\"datac1\")'")
             set_tests_properties(test_OFDM_modem_datac1_octave PROPERTIES PASS_REGULAR_EXPRESSION "BER..: 0.0000")
    # C Tx, C Rx, uncoded       
    add_test(NAME test_OFDM_modem_datac1
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ofdm_mod   --mode datac1 --in  /dev/zero --testframes 10 --verbose 1 |
                            ./ofdm_demod --mode datac1 --out /dev/null --testframes   --verbose 1")
    # C Tx, C Rx, coded       
    add_test(NAME test_OFDM_modem_datac1_ldpc
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ofdm_mod   --mode datac1 --in  /dev/zero --testframes 10 --ldpc --verbose 1 | 
                            ./ofdm_demod --mode datac1 --out /dev/null  --testframes    --ldpc --verbose 1")
           
    # -------------------------------------------------------------------------
    #                                LDPC
    # -------------------------------------------------------------------------

    # tests ldpc_enc/ldpc_noise/ldpc_dec
    add_test(NAME test_ldpc_enc_dec
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ldpc_enc /dev/zero - --sd --code HRA_112_112 --testframes 200 | 
                            ./ldpc_noise - - 0.5 | 
                            ./ldpc_dec - /dev/null --code HRA_112_112 --sd --testframes"
             )

    add_test(NAME test_ldpc_enc_dec_HRAb_396_504
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ldpc_enc /dev/zero - --sd --code HRAb_396_504 --testframes 200 | 
                            ./ldpc_noise - - -2.0 | 
                            ./ldpc_dec - /dev/null --code HRAb_396_504 --sd --testframes"
             )

    add_test(NAME test_ldpc_enc_dec_H_256_768_22
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ldpc_enc /dev/zero - --sd --code H_256_768_22 --testframes 200 | 
                            ./ldpc_noise - - 3.0 | 
                            ./ldpc_dec - /dev/null --code H_256_768_22 --sd --testframes"
             )

    add_test(NAME test_ldpc_enc_dec_H_256_512_4
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ldpc_enc /dev/zero - --sd --code H_256_512_4 --testframes 200 | 
                            ./ldpc_noise - - 0.5 | 
                            ./ldpc_dec - /dev/null --code H_256_512_4 --sd --testframes"
             )

    add_test(NAME test_ldpc_enc_dec_HRAa_1536_512
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ldpc_enc /dev/zero - --sd --code HRAa_1536_512 --testframes 200 | 
                            ./ldpc_noise - - -2 | 
                            ./ldpc_dec - /dev/null --code HRAa_1536_512 --sd --testframes"
             )

    add_test(NAME test_ldpc_enc_dec_H_128_256_5
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ldpc_enc /dev/zero - --sd --code H_128_256_5 --testframes 200 | 
                            ./ldpc_noise - - 0.5 | 
                            ./ldpc_dec - /dev/null --code H_128_256_5 --sd --testframes"
             )

    add_test(NAME test_freedv_api_1600
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_tx 1600 ../../raw/ve9qrp_10s.raw - | ./freedv_rx 1600 - /dev/null")
             set_tests_properties(test_freedv_api_1600 PROPERTIES PASS_REGULAR_EXPRESSION "frames decoded: 503")
           
    add_test(NAME test_freedv_api_700C
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_tx 700C ../../raw/ve9qrp_10s.raw - | ./freedv_rx 700C - /dev/null")
             set_tests_properties(test_freedv_api_700C PROPERTIES PASS_REGULAR_EXPRESSION "frames decoded: 125")
           
    add_test(NAME test_freedv_api_700D_backwards_compatability
             COMMAND sh -c "$<TARGET_FILE:freedv_rx> 700D ${CMAKE_CURRENT_SOURCE_DIR}/raw/testframes_700d.raw /dev/null --testframes --discard"
             )

    # no random speech output due to trial sync when listening to noise
    add_test(NAME test_freedv_api_700D_burble
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_tx 700D ../../raw/ve9qrp.raw - | 
                            ./cohpsk_ch - - -10 --Fs 8000 | 
                            ./freedv_rx 700D - /dev/null --squelch -2")
             set_tests_properties(test_freedv_api_700D_burble PROPERTIES PASS_REGULAR_EXPRESSION "frames decoded: 746  output speech samples: 0")
           
    add_test(NAME test_freedv_api_700D_AWGN_BER
             COMMAND sh -c "dd bs=2560 count=120 if=/dev/zero | $<TARGET_FILE:freedv_tx> 700D - - --testframes | $<TARGET_FILE:cohpsk_ch> - - -20 --Fs 8000 -f -10 | $<TARGET_FILE:freedv_rx> 700D - /dev/null --testframes --discard"
             )

    add_test(NAME test_freedv_api_700D_AWGN_BER_USECOMPLEX
             COMMAND sh -c "dd bs=2560 count=120 if=/dev/zero | $<TARGET_FILE:freedv_tx> 700D - - --testframes | $<TARGET_FILE:cohpsk_ch> - - -20 --Fs 8000 -f -10 | $<TARGET_FILE:freedv_rx> 700D - /dev/null --testframes --discard --usecomplex"
             )

if(LPCNET)

    add_test(NAME test_freedv_api_2020_to_ofdm_demod
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                               ./freedv_tx 2020 ../../wav/wia_16kHz.wav - --testframes | 
                               ./ofdm_demod --mode 2020 --verbose 1 --ldpc -p 312 --testframes > /dev/null"
             )

    add_test(NAME test_freedv_api_2020_from_ofdm_mod
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ofdm_mod --in /dev/zero --mode 2020 --verbose 1 --ldpc -p 312 --testframes 10 | 
                            ./freedv_rx 2020 - /dev/null --testframes"
             )

    add_test(NAME test_freedv_api_2020_awgn
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            dd bs=32000 count=10 if=/dev/zero |
                            ./freedv_tx 2020 - - --testframes | 
                            ./cohpsk_ch - - -24 --Fs 8000 |
                            ./freedv_rx 2020 - /dev/null --testframes"
             )
endif()
           
     add_test(NAME test_freedv_api_2400A
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_tx 2400A ../../raw/ve9qrp_10s.raw - | ./freedv_rx 2400A - /dev/null")
             set_tests_properties(test_freedv_api_2400A PROPERTIES PASS_REGULAR_EXPRESSION "frames decoded: 250")
     add_test(NAME test_freedv_api_2400B
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_tx 2400B ../../raw/ve9qrp_10s.raw - | ./freedv_rx 2400B - /dev/null")
             set_tests_properties(test_freedv_api_2400B PROPERTIES PASS_REGULAR_EXPRESSION "frames decoded: 250")
     add_test(NAME test_freedv_api_800XA
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_tx 800XA ../../raw/ve9qrp_10s.raw - | ./freedv_rx 800XA - /dev/null")
             set_tests_properties(test_freedv_api_800XA PROPERTIES PASS_REGULAR_EXPRESSION "frames decoded: 125")

     add_test(NAME test_freedv_api_rawdata_800XA
             COMMAND sh -c "./tfreedv_800XA_rawdata"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unittest
             )

     add_test(NAME test_freedv_api_rawdata_2400A
             COMMAND sh -c "./tfreedv_2400A_rawdata"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unittest
             )

     add_test(NAME test_freedv_api_rawdata_2400B
             COMMAND sh -c "./tfreedv_2400B_rawdata"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unittest
             )

if(UNIX) # Uses pthreads
     add_test(NAME test_fifo
             COMMAND $<TARGET_FILE:tfifo>
             )
endif()

if (NOT APPLE)
     add_test(NAME test_memory_leak_FreeDV_1600_tx
             COMMAND sh -c " valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./freedv_tx 1600 ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw /dev/null"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
             set_tests_properties(test_memory_leak_FreeDV_1600_tx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")

     add_test(NAME test_memory_leak_FreeDV_1600_rx
             COMMAND sh -c "./freedv_tx 1600 ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw t.raw; \
                            valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./freedv_rx 1600 t.raw /dev/null"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
             set_tests_properties(test_memory_leak_FreeDV_1600_rx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")

     add_test(NAME test_memory_leak_FreeDV_700D_tx
             COMMAND sh -c " valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./freedv_tx 700D ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw /dev/null"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
             set_tests_properties(test_memory_leak_FreeDV_700D_tx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")

     add_test(NAME test_memory_leak_FreeDV_700D_rx
             COMMAND sh -c "./freedv_tx 700D ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw t.raw; \
                            valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./freedv_rx 700D t.raw /dev/null"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
             set_tests_properties(test_memory_leak_FreeDV_700D_rx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")

     add_test(NAME test_memory_leak_FreeDV_700C_tx
             COMMAND sh -c " valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./freedv_tx 700C ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw /dev/null"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
             set_tests_properties(test_memory_leak_FreeDV_700C_tx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")

     add_test(NAME test_memory_leak_FreeDV_700C_rx
             COMMAND sh -c "cd  ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_tx 700C ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw t.raw; \
                            valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./freedv_rx 700C t.raw /dev/null"
             )
             set_tests_properties(test_memory_leak_FreeDV_700C_rx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")

     add_test(NAME test_memory_leak_FreeDV_FSK_LDPC_tx
             COMMAND sh -c "cd  ${CMAKE_CURRENT_BINARY_DIR}/src;
                            valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes \
                            ./freedv_data_raw_tx --testframes 10 FSK_LDPC /dev/zero /dev/null")
             set_tests_properties(test_memory_leak_FreeDV_FSK_LDPC_tx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")

if(LPCNET)
     add_test(NAME test_memory_leak_FreeDV_2020_tx
             COMMAND sh -c "cd  ${CMAKE_CURRENT_BINARY_DIR}/src;
                            valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./freedv_tx 2020 ../../wav/wia_16kHz.wav /dev/null"
             )
             set_tests_properties(test_memory_leak_FreeDV_2020_tx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")

     add_test(NAME test_memory_leak_FreeDV_2020_rx
             COMMAND sh -c "cd  ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_tx 2020 ../../wav/wia_16kHz.wav t.raw; \
                            valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./freedv_rx 2020 t.raw /dev/null"
             )
             set_tests_properties(test_memory_leak_FreeDV_2020_rx PROPERTIES PASS_REGULAR_EXPRESSION "ERROR SUMMARY: 0 errors")
endif(LPCNET)
endif(NOT APPLE)

     add_test(NAME test_codec2_mode_dot_c2
             COMMAND sh -c "./c2enc 700C ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw hts1a.c2 && ./c2dec 1600 hts1a.c2 /dev/null"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
             set_tests_properties(test_codec2_mode_dot_c2 PROPERTIES PASS_REGULAR_EXPRESSION "mode 8")

     add_test(NAME test_codec2_mode_3200
             COMMAND sh -c "./c2enc 3200 ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw - | ./c2dec 3200 - - | sox -t .s16 -r 8000 - hts1a_3200.wav"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
     add_test(NAME test_codec2_mode_2400
             COMMAND sh -c "./c2enc 2400 ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw - | ./c2dec 2400 - - | sox -t .s16 -r 8000 - hts1a_2400.wav"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
     add_test(NAME test_codec2_mode_1400
             COMMAND sh -c "./c2enc 1400 ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw - | ./c2dec 1400 - - | sox -t .s16 -r 8000 - hts1a_1400.wav"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
     add_test(NAME test_codec2_mode_1300
             COMMAND sh -c "./c2enc 1300 ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw - | ./c2dec 1300 - - | sox -t .s16 -r 8000 - hts1a_1300.wav"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
     add_test(NAME test_codec2_mode_1200
             COMMAND sh -c "./c2enc 1200 ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw - | ./c2dec 1200 - - | sox -t .s16 -r 8000 - hts1a_1200.wav"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
     add_test(NAME test_codec2_mode_700C
             COMMAND sh -c "./c2enc 700C ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw - | ./c2dec 700C - - | sox -t .s16 -r 8000 - hts1a_700C.wav"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
     add_test(NAME test_codec2_mode_450
             COMMAND sh -c "./c2enc 450 ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw - | ./c2dec 450 - - | sox -t .s16 -r 8000 - hts1a_450.wav"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )
     add_test(NAME test_codec2_mode_450PWB
             COMMAND sh -c "./c2enc 450PWB ${CMAKE_CURRENT_SOURCE_DIR}/raw/hts1a.raw - | ./c2dec 450PWB - - | sox -t .s16 -r 16000 - hts1a_450PWB.wav"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
             )

     configure_file(unittest/est_n0.sh.in unittest/est_n0.sh)
     add_test(NAME test_est_n0
             COMMAND ${CMAKE_CURRENT_BINARY_DIR}/unittest/est_n0.sh
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/unittest
             )
     add_test(NAME test_vq_mbest
             COMMAND sh -c "./tvq_mbest; \
                            cat target.f32 | ../misc/vq_mbest -k 2 -q vq1.f32,vq2.f32 --mbest 2 -v > out.f32; \
                            diff target.f32 out.f32"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unittest
             )
     configure_file(unittest/test_700c_eq.sh.in unittest/test_700c_eq.sh)
     add_test(NAME test_700c_eq
             COMMAND ${CMAKE_CURRENT_BINARY_DIR}/unittest/test_700c_eq.sh
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/unittest
             )

     # -------------------------------------------------------------------------
     #                                FSK Modem      
     # -------------------------------------------------------------------------
           
     # Octave FSK Modem, to make sure we don't break reference simulation
     add_test(NAME test_fsk_lib
             COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/octave; DISPLAY=\"\" octave --no-gui -qf fsk_lib_demo.m")
             set_tests_properties(test_fsk_lib PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

     add_test(NAME test_fsk_modem_octave_port
             COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/octave; 
                            PATH_TO_TFSK=${CMAKE_CURRENT_BINARY_DIR}/unittest/tfsk octave --no-gui -qf tfsk.m")
             set_tests_properties(test_fsk_modem_octave_port PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

     add_test(NAME test_fsk_modem_mod_demod
             COMMAND sh -c "$<TARGET_FILE:fsk_get_test_bits> - 10000 | 
                            $<TARGET_FILE:fsk_mod> 2 8000 100 1200 100 - - | 
                            $<TARGET_FILE:fsk_demod> -l 2 8000 100 - - | 
                            $<TARGET_FILE:fsk_put_test_bits> -p 99 -q -"
             )
 
     # 2FSK modem at Eb/No = 9dB, SNR = Eb/No+10log10(Rb/B) = 9 + 10*log10(100/3000) = -5.7dB
     # Ideal BER = 0.0094, set thresh 50% higher
     add_test(NAME test_fsk_2fsk_ber
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src; 
                            ./fsk_get_test_bits - 10000 | ./fsk_mod 2 8000 100 1000 100 - - | 
                            ./cohpsk_ch - - -26 --Fs 8000 | 
                            ./fsk_demod 2 8000 100 - - | ./fsk_put_test_bits -b 0.015 -q - ")
     # 4FSK modem at Eb/No = 6dB, SNR = Eb/No+10log10(Rb/B) = 6 + 10*log10(2*100/3000) = -5.7dB
     # Ideal BER = 0.016, set thresh 50% higher
     add_test(NAME test_fsk_4fsk_ber
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src; 
                            ./fsk_get_test_bits - 10000 | ./fsk_mod 4 8000 100 1000 100 - - | 
                            ./cohpsk_ch - - -26 --Fs 8000 | 
                            ./fsk_demod 4 8000 100 - - | ./fsk_put_test_bits -b 0.025 - ")
     # shift FSK signal to -ve frequencies, note No 3dB higher as noise is single sided
     add_test(NAME test_fsk_4fsk_ber_negative_freq
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src; 
                            ./fsk_get_test_bits - 10000 | ./fsk_mod 4 8000 100 1000 200 - - | 
                            ./cohpsk_ch - - -23 --Fs 8000 --ssbfilt 0 --complexout -f -4000 | 
                            ./fsk_demod -c -p 8 4 8000 100 - - | 
                            ./fsk_put_test_bits -b 0.025 -q - ")
     # Low SNR 4FSK uncoded PER/BER test:
     # 4FSK modem at Eb/No = 2dB, SNR = Eb/No+10log10(Rb/B) = 6 + 10*log10(2*100/3000) = -15.7dB
     # Theoretical BER is 0.14.  As it's single sided noise cohpsk_ch reports 3dB less (-18.7dB)
     # Pass condition is 10% PER
     add_test(NAME test_fsk_4fsk_lockdown
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src; 
                            bits=512; tx_packets=20; rx_packets=18; tx_tone_sep=270; Rs=25; 
                            ./fsk_get_test_bits - $(($bits*$tx_packets)) $bits | 
                            ./fsk_mod 4 8000 $Rs 1000 $tx_tone_sep - - | 
                            ./cohpsk_ch - - -13 --Fs 8000 --ssbfilt 0 -f -3000 --complexout | 
                            ./fsk_demod -c -p 8 --mask $tx_tone_sep -t1 --nsym 100 4 8000 $Rs - - 2>stats.txt | 
                            ./fsk_put_test_bits -t 0.25 -b 0.20 -p $rx_packets -f $bits -q -")

     # Octave 4FSK LLR reference simulation - make sure this keeps working
     add_test(NAME test_fsk_lib_4fsk_ldpc
             COMMAND sh -c "cd ${CMAKE_CURRENT_SOURCE_DIR}/octave; DISPLAY=\"\" octave --no-gui -qf fsk_lib_ldpc_demo.m")
             set_tests_properties(test_fsk_lib_4fsk_ldpc PROPERTIES PASS_REGULAR_EXPRESSION "PASS")

     # Command line Unique Word (UW) framer in hard decision mode
     add_test(NAME test_fsk_framer
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src; 
                            ./fsk_get_test_bits - 300 | 
                            ./framer - - 100 51 | 
                            ./deframer - - 100 51 --hard | 
                            ./fsk_put_test_bits -")
             set_tests_properties(test_fsk_framer PROPERTIES PASS_REGULAR_EXPRESSION "PASS")
            
     # Command line Unique Word (UW) framer with LLRs and LDPC (no noise)
     add_test(NAME test_fsk_framer_ldpc
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src; 
                            ./ldpc_enc /dev/zero - --code HRA_112_112 --testframes 10 | ./framer - - 224 51 | 
                            ./tollr | ./deframer - - 224 51 | ./ldpc_dec - /dev/null --code HRA_112_112 --testframes")

     # mFSK soft decision rx_filter to LLR mapping
     add_test(NAME test_fsk_llr
             COMMAND sh -c "${CMAKE_CURRENT_BINARY_DIR}/unittest/tfsk_llr")

     # 4FSK LDPC modem with framer at Rs=100 (uncoded Rb=200), rate 0.8 code
     # SNR = Eb/No + 10*log10(Rb/B) = 5 + 10*log10(200/3000) = -6.7dB
     # Coded Ebc/No = Eb/No - 10*log1010(0.8) = 5 - 10*log10(0.8) = 6.0dB
     # (calculation ignores small UW overhead).  See also test_freedv_fsk_ldpc below
     # which is the same thing bundled up into a FreeDV "mode"
     add_test(NAME test_fsk_4fsk_ldpc
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./ldpc_enc /dev/zero - --code HRAb_396_504 --testframes 200 | 
                            ./framer - - 504 5186 |
                            ./fsk_mod 4 8000 100 1000 100 - - | 
                            ./cohpsk_ch - - -25 --Fs 8000  | 
                            ./fsk_demod -s 4 8000 100 - - | 
                            ./deframer - - 504 5186  | 
                            ./ldpc_dec - /dev/null --code HRAb_396_504 --testframes")

     # 800XA framer test
     add_test(NAME test_fsk_vhf_framer
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src; 
                            ./c2enc 700C ../../raw/ve9qrp_10s.raw - | 
                            ./vhf_frame_c2 B - - | 
                            ./fsk_mod -p 10 4 8000 400 400 400 - - | 
                            ./fsk_demod -p 10 4 8000 400 - - | 
                            ./vhf_deframe_c2 B - /dev/null")
             set_tests_properties(test_fsk_vhf_framer PROPERTIES PASS_REGULAR_EXPRESSION "total_uw_err: 0")

     # VHF Ethernet-style packet system
     add_test(NAME test_freedv_data_channel
             COMMAND sh -c "./tfreedv_data_channel"
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/unittest
             )
           
     # ---------------------------------------------------------
     # FreeDV API raw data
     # ---------------------------------------------------------
     
     # OFDM LDPC
     add_test(NAME test_freedv_data_raw_ofdm
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            head -c 140 </dev/urandom > binaryIn.bin;
                            ./freedv_data_raw_tx 700D binaryIn.bin - | 
                            ./freedv_data_raw_rx 700D - - -v | 
                            diff /dev/stdin binaryIn.bin")

     # FSK LDPC default 100 bit/s 2FSK, enough noise for several % raw BER to give
     # FEC/acquisition a work out, bursts of 1 frame as that stresses acquisition
     add_test(NAME test_freedv_data_raw_fsk_ldpc_100
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_data_raw_tx --testframes 1 --bursts 10 FSK_LDPC /dev/zero - |
                            ./cohpsk_ch - - -5 --Fs 8000 --ssbfilt 0 |
                            ./freedv_data_raw_rx --testframes -v FSK_LDPC - /dev/null")           
             set_tests_properties(test_freedv_data_raw_fsk_ldpc_100 PROPERTIES PASS_REGULAR_EXPRESSION "output_packets: 10")

     # FSK LDPC 1000 bit/s 2FSK, Fs=40kHz, as different configs can upset acquisition
     add_test(NAME test_freedv_data_raw_fsk_ldpc_1k
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_data_raw_tx --Fs 40000 --Rs 1000 --tone1 1000 --shift 1000 --testframes 1 --bursts 10 FSK_LDPC /dev/zero - |
                            ./cohpsk_ch - - -10 --Fs 8000 --ssbfilt 0 |
                            ./freedv_data_raw_rx --testframes -v --Fs 40000 --Rs 1000 FSK_LDPC - /dev/null")
             set_tests_properties(test_freedv_data_raw_fsk_ldpc_1k PROPERTIES PASS_REGULAR_EXPRESSION "output_packets: 10")

     # FSK LDPC 10000 bit/s 2FSK, Fs=100kHz, each of the 10 bursts has 100 frames
     add_test(NAME test_freedv_data_raw_fsk_ldpc_10k
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_data_raw_tx --Fs 100000 --Rs 10000 --tone1 10000 --shift 10000 --testframes 100 --bursts 10 FSK_LDPC /dev/zero - |
                            ./cohpsk_ch - - -15 --Fs 8000 --ssbfilt 0 | 
                            ./freedv_data_raw_rx --testframes -v --Fs 100000 --Rs 10000 FSK_LDPC - /dev/null")
             set_tests_properties(test_freedv_data_raw_fsk_ldpc_10k PROPERTIES PASS_REGULAR_EXPRESSION "output_packets: 1000")

     # FSK LDPC Rs=1000 bit/s (Rb=2000) 4FSK, Fs=40kHz, this needs --mask and longer preamble to work, due to lower Es/No, about 2dB over 2FSK
     add_test(NAME test_freedv_data_raw_fsk_ldpc_2k
             COMMAND sh -c "cd ${CMAKE_CURRENT_BINARY_DIR}/src;
                            ./freedv_data_raw_tx -a 8192 -m 4 --Fs 40000 --Rs 1000 --tone1 1000 --shift 1000 --testframes 1 --bursts 10 FSK_LDPC /dev/zero - |
                             ./cohpsk_ch - - -22 --Fs 8000 --ssbfilt 0 | 
                             ./freedv_data_raw_rx -m 4 --testframes -v --Fs 40000 --Rs 1000 FSK_LDPC --mask 1000 - /dev/null")
             set_tests_properties(test_freedv_data_raw_fsk_ldpc_2k PROPERTIES PASS_REGULAR_EXPRESSION "output_packets: 10")
             
endif(UNITTEST)
