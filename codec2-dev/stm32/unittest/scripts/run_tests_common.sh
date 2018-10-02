# This file must be "sourced" from a parent shell!
#
# run_tests_common.sh
#
# This is a collection of common variable settings for stm32 unit tests.
#
# The variable $SCRIPTS must be set when this is called.

if [ -z ${SCRIPTS+x} ]; then 
    echo "Error, run_tests_common.sh requires that \$SCRIPTS be set!"
    exit 1
    fi

#######################################
# Set default directories based on the parent of the SCRIPTS variable.
set -a 

#UNITTEST_BASE - Location of STM32 Unittests and files
UNITTEST_BASE="$( cd "$( dirname "${SCRIPTS}" )" >/dev/null && pwd )"

# UNITTEST_BIN - Location of STM32 unittest binaries
UNITTEST_BIN="${UNITTEST_BASE}/src"

# STM32_BASE - Base directory of Codec2
STM32_BASE="$( cd "$( dirname "${UNITTEST_BASE}" )" >/dev/null && pwd )"

# CODEC2_BASE - Base directory of Codec2
CODEC2_BASE="$( cd "$( dirname "${STM32_BASE}" )" >/dev/null && pwd )"

# CODEC2_BIN - Location of x86 utiliy programs for Codec2
CODEC2_BIN="${CODEC2_BASE}/build_linux/src"

# CODEC2_UTST - Location of x86 utiliy programs for Codec2 unittest
CODEC2_UTST="${CODEC2_BASE}/build_linux/unittest"

set +a 

#######################################
# Add directories to PATH
export PATH=${PATH}:${SCRIPTS}:${CODEC2_BIN}:${CODEC2_UTST}


#######################################
# Parse command line options
# Options (starting with "--") are stored in $ARGS.
# Non-options are taken as the test name (last one sticks).
declare -A ARGS
for arg in "$@"; do
    if [[ ${arg} == --* ]] ; then ARGS[${arg}]=true
    else TEST=${arg}
    fi
    done

#######################################
# A function for setup

setup_common () {

    if [ ${ARGS[--clean]+_} ] ; then
        if [ -d "${1}" ] ; then rm -rf "${1}"; fi
        fi

    # Make run directory if needed
    if [ ! -d "${1}" ] ; then mkdir -p "${1}"; fi

    }
