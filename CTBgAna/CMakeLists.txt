CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
include_directories(src)

set(SRC main.cpp src/cosmicTemplatesBgAna.cpp ../CTCosPatTrain/cosmicTemplatesEfficiency.cpp)
add_executable(Mu3eCTBg2Ana ${SRC})
target_link_libraries(Mu3eCTBg2Ana patterns rootTreeFiles)

set(SRC main_manybgdata.cpp src/cosmicTemplatesBgAna.cpp ../CTCosPatTrain/cosmicTemplatesEfficiency.cpp)
add_executable(Mu3eCTBg2MultiAna ${SRC})
target_link_libraries(Mu3eCTBg2MultiAna patterns rootTreeFiles)

set(SRC main_bgeffplot.cpp src/cosmicTemplatesBgPlots.cpp)
add_executable(Mu3eCTBg3Plots ${SRC})
target_link_libraries(Mu3eCTBg3Plots rootTreeFiles)

set(SRC main_CTBkgAnaDatapoints.cpp src/cosmicTemplatesBgAna.cpp ../CTCosPatTrain/cosmicTemplatesEfficiency.cpp)
add_executable(Mu3eCTBg2AnaDatapoints ${SRC})
target_link_libraries(Mu3eCTBg2AnaDatapoints patterns rootTreeFiles)