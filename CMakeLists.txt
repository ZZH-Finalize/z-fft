cmake_minimum_required(VERSION 3.20)
project(z_fft LANGUAGES C)

# 如果更上层的构建系统指定了基本数据类型, 则使用更上层的, 否则默认使用float
if(NOT ${ZFFT_BASIC_TYPE} STREQUAL "")
    message("Using external type:" ${ZFFT_BASIC_TYPE})
    add_definitions(-DBASIC_TYPE=${ZFFT_BASIC_TYPE})
else()
    add_definitions(-DBASIC_TYPE=float)
endif()

include_directories(.)

# 本项目核心, 这里编译了一个共享库, 这实际上是本项目的实际输出, libfft.dll或libfft.so
add_library(fft SHARED
    fft.c
)

# 这里并不是真的需要一个可执行程序, 只是编译了一个测试案例而已, 主要还是提供上面那个库
add_executable(${PROJECT_NAME}
    main.c
)

target_link_libraries(${PROJECT_NAME} fft)
