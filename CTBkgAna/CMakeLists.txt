CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
include_directories(src)

set(SRC main_CTBkgAnaMany.cpp src/cosmicTemplatesBgAna.cpp ../CTCosPatTrain/src/getCosmicSIDtracks.cpp)
add_executable(EvalBkgDBParamCombinations ${SRC})
target_link_libraries(EvalBkgDBParamCombinations patterns rootTreeFiles)

set(SRC main_CTBkgAnaDatapoints.cpp src/cosmicTemplatesBgAna.cpp ../CTCosPatTrain/src/getCosmicSIDtracks.cpp)
add_executable(EvalBkgDBMulti ${SRC})
target_link_libraries(EvalBkgDBMulti patterns rootTreeFiles)