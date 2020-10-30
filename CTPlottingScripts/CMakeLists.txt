CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
include_directories(src)

set(SRC2 main_bgAnaROC.cpp src/bgAnaPlots_ROC.cpp)
add_executable(PlotsBgAnaROC ${SRC2})
target_link_libraries(PlotsBgAnaROC rootTreeFiles)

set(SRC3 main_bgAnaEffSPcount.cpp src/bgAnaPlots_EffSPcount.cpp)
add_executable(PlotsBgAnaSPcount ${SRC3})
target_link_libraries(PlotsBgAnaSPcount rootTreeFiles)

set(SRC6 main_bgAnaSPC.cpp src/bgAnaPlots_EffSPC.cpp)
add_executable(PlotsBgAnaSPCdatapoints ${SRC6})
target_link_libraries(PlotsBgAnaSPCdatapoints rootTreeFiles)

set(SRC7 main_bgAnaSPR.cpp src/bgAnaPlots_EffSPR.cpp)
add_executable(PlotsBgAnaSPRdatapoints ${SRC7})
target_link_libraries(PlotsBgAnaSPRdatapoints rootTreeFiles)

set(SRC5 main_bgAnaEffBeam.cpp src/bgAnaPots_EffBeamRate.cpp)
add_executable(PlotsBgAnaEffBeam ${SRC5})
target_link_libraries(PlotsBgAnaEffBeam rootTreeFiles)

set(SRC4 main_CTtrainplots.cpp src/cosmicTemplatesTrainingPlots.cpp)
add_executable(PlotsCosPatTrain ${SRC4})
target_link_libraries(PlotsCosPatTrain rootTreeFiles)

set(SRC5 main_CTtrainTemplSPC.cpp src/ctTrainingPlots_TemplatesSPC.cpp)
add_executable(PlotsCosPatTemplSPC ${SRC5})
target_link_libraries(PlotsCosPatTemplSPC rootTreeFiles)

