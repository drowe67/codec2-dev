#INCLUDE(CMakeForceCompiler)
 
SET(CMAKE_SYSTEM_NAME GNU)
SET(CMAKE_SYSTEM_VERSION 1)
SET(CMAKE_CROSSCOMPILING TRUE)
 
# specify the cross compiler
set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-cpp)
 
#SET(COMMON_FLAGS "-mcpu=cortex-m3 -mthumb -mthumb-interwork -msoft-float -ffunction-sections -fdata-sections -g -fno-common -fmessage-length=0")
#SET(CMAKE_CXX_FLAGS "${COMMON_FLAGS}  -std=gnu++0x")
#SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu99")
#set(CMAKE_EXE_LINKER_FLAGS "-Wl,-gc-sections ")
