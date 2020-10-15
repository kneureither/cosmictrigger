CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
include_directories(src)

#add_library(cospat STATIC processSegsPrototype.cpp SuperPixelHit.cpp PatternEngineSingle.cpp)
#add_library(util SHARED ../util/plots.h ../util/utilityFunctions.h)
#SET_TARGET_PROPERTIES(util PROPERTIES LINKER_LANGUAGE C)
#target_link_libraries(cospat ${ROOT_LIBRARIES})

set(SRC main_CTbuildMany.cpp cosmicTemplatesBuild.cpp cosmicTemplatesTrainingPlots.cpp)
add_executable(Mu3eCosPat1BuildMany ${SRC})
target_link_libraries(Mu3eCosPat1BuildMany patterns rootTreeFiles)

set(SRC2 main_CTbuild.cpp cosmicTemplatesBuild.cpp ../root/TemplateDatabaseFile.cpp)
add_executable(Mu3eCosPatBuild ${SRC2})
target_link_libraries(Mu3eCosPatBuild patterns rootTreeFiles)

set(SRC3 main_CTefficiency.cpp cosmicTemplatesEfficiency.cpp ../root/TemplateDatabaseFile.cpp)
add_executable(Mu3eCosPat3Eff ${SRC3})
target_link_libraries(Mu3eCosPat3Eff patterns rootTreeFiles)

set(SRC4 main_CTplots.cpp cosmicTemplatesTrainingPlots.cpp)
add_executable(Mu3eCosPat2Plots ${SRC4})
target_link_libraries(Mu3eCosPat2Plots rootTreeFiles)

set(SRC5 main_get_cosmic_data.cpp cosmicTemplatesEfficiency.cpp)
add_executable(Mu3eCosPat4LoadCosData ${SRC5})
target_link_libraries(Mu3eCosPat4LoadCosData rootTreeFiles patterns)