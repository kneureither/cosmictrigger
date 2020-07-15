CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

#add_library(cospat STATIC processSegsPrototype.cpp SuperPixelHit.cpp PatternEngineSingle.cpp)
#add_library(util SHARED ../util/plots.h ../util/utilityFunctions.h)
#SET_TARGET_PROPERTIES(util PROPERTIES LINKER_LANGUAGE C)
#target_link_libraries(cospat ${ROOT_LIBRARIES})

add_library(patterns SHARED SuperPixelHit.cpp PatternEngineSingle.cpp TemplateBank.cpp TemplateData.cpp PatternEngine.cpp SPCalculations.cpp)
target_link_libraries(patterns ${ROOT_LIBRARIES})

set(SRC main.cpp processSegsPrototype.cpp)
add_executable(Mu3eCosPatPrototype ${SRC})
target_link_libraries(Mu3eCosPatPrototype patterns)

set(SRC2 main.cpp buildCosmicTrajTemplatesScript.cpp ../root/TemplateDatabaseTree.cpp)
add_executable(Mu3eCosPat ${SRC2})
target_link_libraries(Mu3eCosPat patterns rootTreeFiles)

#target_link_libraries(Mu3eCosPat ${ROOT_LIBRARIES} f2clib karilib)
#add_compile_options("-lf2c" "-DMAIN")