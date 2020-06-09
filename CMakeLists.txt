CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

set(F2CPATH /usr/local/include/f2c.h)
include_directories(${F2CPATH})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

set(KARIPATH "../3rdparty/karimaki/")
set(SRC
        scripts/recacc_main.cpp
        scripts/reconstruction_accuracy.cpp
        ${KARIPATH}karimaki_hit.c
        ${KARIPATH}fitszw.c
        ${KARIPATH}cirpar.c
        ${KARIPATH}cirparw.c
        ${KARIPATH}wrappers.c
        )


add_executable(Mu3eRecAcc ${SRC})

add_compile_options("-lf2c" "-DMAIN")

add_subdirectory(scripts)
add_subdirectory(tests)