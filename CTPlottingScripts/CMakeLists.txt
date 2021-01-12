CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
include_directories(src)

set(SRC2 main_bgAnaPlots_legacyROC.cpp src/bgAnaPlots_legacyROC.cpp)
add_executable(PlotsBgAnaLegacyROC ${SRC2})
target_link_libraries(PlotsBgAnaLegacyROC rootTreeFiles)

set(SRC3 main_bgAnaPlots_legacySPC.cpp src/bgAnaPlots_legacySPC.cpp)
add_executable(PlotsBgAnaLegacySPC ${SRC3})
target_link_libraries(PlotsBgAnaLegacySPC rootTreeFiles)

set(SRC6 main_bgAnaPlots_SPC.cpp src/bgAnaPlots_EffSPC.cpp)
add_executable(PlotsBgAnaSPC ${SRC6})
target_link_libraries(PlotsBgAnaSPC rootTreeFiles)

set(SRC7 main_bgAnaPlots_SPR.cpp src/bgAnaPlots_EffSPR.cpp)
add_executable(PlotsBgAnaSPR ${SRC7})
target_link_libraries(PlotsBgAnaSPR rootTreeFiles)

set(SRC7 main_bgAnaPlots_ROC.cpp src/bgAnaPlots_EffROC.cpp)
add_executable(PlotsBgAnaROC ${SRC7})
target_link_libraries(PlotsBgAnaROC rootTreeFiles)

set(SRC5 main_bgAnaPlots_EffBeamRate.cpp src/bgAnaPlots_EffBeamRate.cpp)
add_executable(PlotsBgAnaEffBeam ${SRC5})
target_link_libraries(PlotsBgAnaEffBeam rootTreeFiles)

set(SRC4 main_CTtrainplots.cpp src/ctTrainingPlots.cpp)
add_executable(PlotsDBTraining ${SRC4})
target_link_libraries(PlotsDBTraining rootTreeFiles)

set(SRC5 main_CTtrainTemplSPC.cpp src/ctTrainingPlots_TemplatesSPC.cpp)
add_executable(PlotsBuildDBTemplSPC ${SRC5})
target_link_libraries(PlotsBuildDBTemplSPC rootTreeFiles)

set(SRC5 PlotPEconfig.cpp)
add_executable(PlotsPEconfig ${SRC5})
target_link_libraries(PlotsPEconfig rootTreeFiles)

