### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# CMake configuration
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
cmake_minimum_required(VERSION 3.16)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Project configuration
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
project(experiments VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_VERBOSE_MAKEFILE ON)

# Debug build: Activate GLIBC assertions
add_compile_definitions($<$<CONFIG:DEBUG>:_GLIBCXX_ASSERTIONS>)

add_executable(experiments
        src/utils.hpp
        src/main.cpp
        src/tseries/readers/tsreader/tsreader.hpp
        src/tseries/readers/tsreader/tsreader.cpp
        src/tseries/tseries.hpp
        src/distances/distances.hpp
        src/distances/dtw/wdtw.hpp
        src/distances/dtw/cdtw.hpp
        src/distances/dtw/pruneddtw.hpp
        src/distances/dtw/dtw.hpp
        src/distances/elementwise/elementwise.hpp
        src/distances/lcss/lcss.hpp
        src/distances/twe/twe.hpp
        src/distances/erp/erp.hpp
        src/distances/msm/msm.hpp
        src/distances/dtw/lowerbounds/envelopes.hpp
        src/distances/dtw/lowerbounds/lb_keogh.hpp
        src/distances/dtw/lowerbounds/lb_webb.hpp
        )

