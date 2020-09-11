CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
include_directories(src)

set(SRC main.cpp src/cosmicTemplatesBgEval.cpp)
add_executable(Mu3eCTBg2Ana ${SRC})
target_link_libraries(Mu3eCTBg2Ana patterns rootTreeFiles)

set(SRC main_bgeffplot.cpp src/cosmicTemplatesBgPlots.cpp)
add_executable(Mu3eCTBg3Plots ${SRC})
target_link_libraries(Mu3eCTBg3Plots rootTreeFiles)

set(SRC main_bg_combinatorics.cpp src/bgDataCombinatorics.cpp ../root/BackgroundDataFile.cpp)
add_executable(Mu3eCTBg1ComputeCmbntcs ${SRC})
target_link_libraries(Mu3eCTBg1ComputeCmbntcs rootTreeFiles)