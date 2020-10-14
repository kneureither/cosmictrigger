CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

set(F2CPATH "/usr/local/include")
set(KARIPATH "${CMAKE_CURRENT_SOURCE_DIR}/karimaki/")

set(F2CLIB /usr/local/lib/libf2c.a)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/util)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/karimaki)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/root)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Mu3eCosPat)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Mu3eCosPat/include) #this is for the core modules

include_directories(${F2CPATH}) #findet f2c.h

add_library(f2clib SHARED ${F2CLIB})
SET_TARGET_PROPERTIES(f2clib PROPERTIES LINKER_LANGUAGE C)

add_library(karilib SHARED ${KARIPATH}karimakiHelixfit.cpp ${KARIPATH}fitszw.cpp ${KARIPATH}cirpar.c ${KARIPATH}cirparw.c ${KARIPATH}wrappers.c)
SET_TARGET_PROPERTIES(karilib PROPERTIES LINKER_LANGUAGE CXX)


add_subdirectory(Mu3eRecAcc)
add_subdirectory(Mu3eCosPat)
add_subdirectory(Mu3ePreSlimData)
add_subdirectory(Mu3eCTBgAna)
add_subdirectory(tests)
add_subdirectory(root)
add_subdirectory(playground)