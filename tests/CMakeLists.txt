CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(mu3eRecAcc C CXX)
set(CMAKE_CXX_STANDARD 14)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

set(KARIPATH "../3rdparty/karimaki/")
set(SRC
        tests.cpp
        )


add_executable(Tests ${SRC})

add_compile_options("-DMAIN")