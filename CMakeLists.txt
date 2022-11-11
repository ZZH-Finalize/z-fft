cmake_minimum_required(VERSION 3.20)
project(z_fft LANGUAGES C)

add_definitions(-DBASIC_TYPE=double)

add_library(fft SHARED
    fft.c
)

add_executable(${PROJECT_NAME}
    main.c
)

target_link_libraries(${PROJECT_NAME} fft)
