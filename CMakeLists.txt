CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/util)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/karimaki)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tests)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/root)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CTCosPatTrain)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CTBkgAna/inc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CTCoreModules/inc)

# Change this to the path of your f2c.h file from the f2c library
set(F2CPATH "/usr/local/include")
set(F2CLIB /usr/local/lib/libf2c.a)
include_directories(${F2CPATH}) #finds f2c.h

add_library(f2clib SHARED ${F2CLIB})
SET_TARGET_PROPERTIES(f2clib PROPERTIES LINKER_LANGUAGE C)

set(KARIPATH "${CMAKE_CURRENT_SOURCE_DIR}/karimaki/")
add_library(karilib SHARED ${KARIPATH}karimakiHelixfit.cpp ${KARIPATH}fitszw.cpp ${KARIPATH}cirpar.c ${KARIPATH}cirparw.c ${KARIPATH}wrappers.c)
SET_TARGET_PROPERTIES(karilib PROPERTIES LINKER_LANGUAGE CXX)


add_subdirectory(RecAccCheck)
add_subdirectory(CTCosPatTrain)
add_subdirectory(CTPreSlimData)
add_subdirectory(CTBkgAna)
add_subdirectory(CTPlottingScripts)
add_subdirectory(tests)
add_subdirectory(root)
