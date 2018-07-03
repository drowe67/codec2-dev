#INCLUDE(CMakeForceCompiler)
 
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)
#set(CMAKE_SYSTEM_VERSION 1)
#set(CMAKE_CROSSCOMPILING TRUE)
 
# specify the cross compiler
set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-cpp)

set(CMAKE_EXE_LINKER_FLAGS "-T${CMAKE_SOURCE_DIR}/stm32_flash.ld")
set(CMAKE_EXECUTABLE_SUFFIX_C ".elf")
set(CMAKE_EXECUTABLE_SUFFIX_CXX ".elf")
set(CMAKE_EXECUTABLE_SUFFIX_ASM ".elf")

# Macro for elf->bin
macro(elf2bin target)
    add_custom_command(TARGET ${target}
    POST_BUILD COMMAND ${CMAKE_OBJCOPY} -O binary ${target}.elf ${target}.bin
    COMMENT "Creating binary for ${target}")
endmacro()
