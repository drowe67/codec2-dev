# This file must be "sourced" from a parent shell!
#
# setup.sh
#
# This is a collection of common variable settings for manually running
# stm32 unit tests.
#
# This assumes it is called from the "stm32/unittests" directory!!!

SCRIPTS="${PWD}/scripts"

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

# CODEC2_UTST - Location of codec2 unittest (and its scripts)
CODEC2_UTST="${CODEC2_BASE}/unittest"

# CODEC2_UTST_BIN - Location of x86 utiliy programs for Codec2 unittest
CODEC2_UTST_BIN="${CODEC2_BASE}/build_linux/unittest"

# CODEC2_SCRIPT - Location of Codec2 scripts
CODEC2_SCRIPT="${CODEC2_BASE}/script"

set +a 

#######################################
# Add directories to PATH(s)
export PATH=${PATH}:${SCRIPTS}
export PATH=${PATH}:${CODEC2_BIN}:${CODEC2_UTST}:${CODEC2_UTST_BIN}:${CODEC2_SCRIPT}
export LD_LIBRARY_PATH=${CODEC2_BIN}:${LD_LIBRARY_PATH}
