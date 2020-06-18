CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

set(F2CPATH "/usr/local/include")
set(KARIPATH "../3rdparty/karimaki/")

set(F2CLIB /usr/local/lib/libf2c.a)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/util)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/karimaki)
include_directories(${F2CPATH}) #findet f2c.h

#kann weg
#add_library(f2c SHARED ${F2CLIB})
#SET_TARGET_PROPERTIES(f2c PROPERTIES LINKER_LANGUAGE C)
#target_link_libraries(f2c ${ROOT_LIBRARIES})


#include(${ROOT_USE_FILE})
#include_directories(${ROOT_INCLUDE_DIR})
#link_directories(${ROOT_LIBRARY_DIR})

add_library(fits ${KARIPATH}cirpar.c ${KARIPATH}cirparw.c)
target_link_libraries(fits ${ROOT_LIBRARIES})

target_link_libraries(fits ${ROOT_LIBRARIES})
#SET_TARGET_PROPERTIES(root PROPERTIES LINKER_LANGUAGE CXX)


set(SRC
        Mu3eRecAcc/recacc_main.cpp
        Mu3eRecAcc/reconstruction_accuracy.cpp
        ${KARIPATH}karimaki_hit.cpp
        ${KARIPATH}fitszw.cpp
        ${KARIPATH}cirpar.c
        ${KARIPATH}cirparw.c
        ${KARIPATH}wrappers.c
        )



add_executable(Mu3eRecAcc ${SRC})
#target_link_libraries(Mu3eRecAcc f2c)
target_link_libraries(Mu3eRecAcc fits)
#target_link_libraries(Mu3eRecAcc -L/usr/local/lib/libf2c.a)
#target_compile_options(Mu3eRecAcc PUBLIC -lf2c -DMAIN)
add_compile_options("-lf2c" "-DMAIN")


add_subdirectory(Mu3eRecAcc)
add_subdirectory(Mu3eCosPat)
add_subdirectory(tests)