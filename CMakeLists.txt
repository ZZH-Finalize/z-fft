cmake_minimum_required(VERSION 3.20)
project(z_fft LANGUAGES C)

if(NOT ${ZFFT_BASIC_TYPE} STREQUAL "")
    message("Using extenal type")
    add_definitions(-DBASIC_TYPE=${ZFFT_BASIC_TYPE})
else()
    add_definitions(-DBASIC_TYPE=double)
endif()

include_directories(.)

add_library(fft SHARED
    fft.c
)

add_executable(${PROJECT_NAME}
    main.c
)

target_link_libraries(${PROJECT_NAME} fft)
