set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)
set(CMAKE_SYSTEM_VERSION 1)
 
# specify the cross compiler
set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-cpp)
set(CMAKE_C_FLAGS_INIT "--specs=nosys.specs")
set(CMAKE_CXX_FLAGS_INIT "--specs=nosys.specs")

# Macro for elf->bin
macro(elf2bin target)
    add_custom_command(TARGET ${target}
    POST_BUILD COMMAND ${CMAKE_OBJCOPY} -O binary ${target}.elf ${target}.bin
    COMMENT "Creating binary for ${target}")
endmacro()
