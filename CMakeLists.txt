CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

set(F2CPATH /usr/local/include)
#set(F2CLIB /usr/local/lib/libf2c.a)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/util)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/karimaki)
include_directories(${F2CPATH})

#add_library(f2c SHARED ${F2CLIB})
#target_link_libraries(f2c ${ROOT_LIBRARIES})


#include(${ROOT_USE_FILE})
#include_directories(${ROOT_INCLUDE_DIR})
#link_directories(${ROOT_LIBRARY_DIR})
#

set(KARIPATH "../3rdparty/karimaki/")
set(SRC
        Mu3eRecAcc/recacc_main.cpp
        Mu3eRecAcc/reconstruction_accuracy.cpp
        ${KARIPATH}karimaki_hit.c
        ${KARIPATH}fitszw.c
        ${KARIPATH}cirpar.c
        ${KARIPATH}cirparw.c
        ${KARIPATH}wrappers.c
        )


add_executable(Mu3eRecAcc ${SRC})
target_compile_options(Mu3eRecAcc PUBLIC -lf2c -DMAIN)
#add_compile_options("-lf2c" "-DMAIN")
target_link_libraries(Mu3eRecAcc -L/usr/local/lib/libf2c.a)

add_subdirectory(Mu3eRecAcc)
add_subdirectory(Mu3eCosPat)
add_subdirectory(tests)