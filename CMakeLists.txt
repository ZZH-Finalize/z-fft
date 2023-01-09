cmake_minimum_required(VERSION 3.20)
project(z_fft LANGUAGES C)
enable_testing()

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

# 性能测试案例, 默认进行16384点变换, 需要耗时几十秒才能跑完, 其中主要是dft算的很慢, fft只需要几毫秒
add_executable(benchmark
    benchmark.c
)
target_link_libraries(benchmark fft)

function(add_testcase TARGET SOURCE_FILES)
    add_executable(${TARGET} ${SOURCE_FILES})
    target_link_libraries(${TARGET} fft)
    add_test(${TARGET}_TestCase ${TARGET}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        CONFIGURATIONS Debug
    )
    # set_tests_properties(${TARGET}_TestCase PROPERTIES PASS_REGULAR_EXPRESSION "${OUTPUT_REGEX}")
endfunction()

add_testcase(demo
    demo.c
)

add_testcase(benchmark_fast
    benchmark_fast.c
)

add_testcase(complex_test
    complex_test.c
)
