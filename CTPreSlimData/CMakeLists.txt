CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

#add_library(cospat STATIC processSegsPrototype.cpp SuperPixelHit.cpp PatternEngine.cpp)
#add_library(util SHARED ../util/plots.h ../util/utilityFunctions.h)
#SET_TARGET_PROPERTIES(util PROPERTIES LINKER_LANGUAGE C)
#target_link_libraries(cospat ${ROOT_LIBRARIES})

set(SRC main.cpp slimSegsDataScript.cpp ../root/SegsTreeRead.cpp ../root/SlimSegsTree.cpp)

add_executable(PreSlimCosData ${SRC})
target_link_libraries(PreSlimCosData ${ROOT_LIBRARIES} f2clib karilib)
add_compile_options("-lf2c" "-DMAIN")

set(SRC PlotSlimSegsTrackVars.cpp ../root/SegsTreeRead.cpp ../root/SlimSegsTree.cpp)

add_executable(PreSlimPlot ${SRC})
target_link_libraries(PreSlimPlot ${ROOT_LIBRARIES} f2clib karilib)
add_compile_options("-lf2c" "-DMAIN")