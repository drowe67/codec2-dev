# This file must be "sourced" from a parent shell!
#
# setup.sh
#
# This is a collection of common variable settings for manually running
# stm32 unit tests.
#
# This assumes it is called from the "stm32/unittests" directory!!!

SCRIPTS="${PWD}/scripts"

# Setup common variables
source $SCRIPTS/run_tests_common.sh

#######################################
# Add directories to PATH(s)
export PATH=${PATH}:${SCRIPTS}
export PATH=${PATH}:${CODEC2_BIN}:${CODEC2_UTST}:${CODEC2_UTST_BIN}:${CODEC2_SCRIPT}
export LD_LIBRARY_PATH=${CODEC2_BIN}:${LD_LIBRARY_PATH}
