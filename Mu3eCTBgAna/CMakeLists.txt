CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
include_directories(src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/root)

set(SRC main.cpp src/cosmicTemplatesBgEval.cpp)
add_executable(Mu3eCTBg2Ana ${SRC})
target_link_libraries(Mu3eCTBg2Ana patterns rootTreeFiles)

set(SRC main_bgeffplot.cpp src/cosmicTemplatesBgPlots.cpp)
add_executable(Mu3eCTBg3Plots ${SRC})
target_link_libraries(Mu3eCTBg3Plots rootTreeFiles)

set(SRC2 main_bgroc.cpp src/cosmicTemplatesBgROC.cpp)
add_executable(Mu3eCTBg4ROC ${SRC2})
target_link_libraries(Mu3eCTBg4ROC rootTreeFiles)

set(SRC2 main_bgeffcount.cpp src/BgEvalEff.cpp)
add_executable(Mu3eCTBg5Eff ${SRC2})
target_link_libraries(Mu3eCTBg5Eff rootTreeFiles)