CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(include)

#add_library(cospat STATIC processSegsPrototype.cpp SuperPixelHit.cpp PatternEngineSingle.cpp)
#add_library(util SHARED ../util/plots.h ../util/utilityFunctions.h)
#SET_TARGET_PROPERTIES(util PROPERTIES LINKER_LANGUAGE C)
#target_link_libraries(cospat ${ROOT_LIBRARIES})

add_library(patterns SHARED SuperPixelHit.cpp PatternEngineSingle.cpp TemplateBank.cpp TemplateData.cpp PatternEngine.cpp SPCalculations.cpp ../root/TemplateDatabase.cpp)
target_link_libraries(patterns ${ROOT_LIBRARIES})

set(SRC main.cpp cosmicTemplatesBuild.cpp makePlots.cpp)
add_executable(Mu3eCosPatSerie ${SRC})
target_link_libraries(Mu3eCosPatSerie patterns rootTreeFiles)

set(SRC2 main_CTbuild.cpp cosmicTemplatesBuild.cpp ../root/TemplateDatabase.cpp)
add_executable(Mu3eCosPatBuild ${SRC2})
target_link_libraries(Mu3eCosPatBuild patterns rootTreeFiles)

set(SRC3 main_CTeval.cpp cosmicTemplatesEval.cpp ../root/TemplateDatabase.cpp)
add_executable(Mu3eCosPatEval ${SRC3})
target_link_libraries(Mu3eCosPatEval patterns rootTreeFiles)

set(SRC4 main_CTplots.cpp makePlots.cpp)
add_executable(Mu3eCosPatPlots ${SRC4})
target_link_libraries(Mu3eCosPatPlots rootTreeFiles)

#target_link_libraries(Mu3eCosPat ${ROOT_LIBRARIES} f2clib karilib)
#add_compile_options("-lf2c" "-DMAIN")