CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
add_library(patterns SHARED ../CTCoreModules/src/PatternEngineSingle.cpp
                            ../CTCoreModules/src/TemplateBank.cpp
                            ../CTCoreModules/src/TemplateData.cpp
                            ../CTCoreModules/src/PatternEngine.cpp
                            ../CTCoreModules/src/SPCalculations.cpp
                            ../root/TemplateDatabaseFile.cpp)

target_link_libraries(patterns ${ROOT_LIBRARIES})

# Builds all template banks (template dbs) of combinations from given parameters
# Plotting script is included, as it is called in the very end. Could be removed.
set(SRC main_CTbuildMany.cpp src/cosmicTemplatesBuild.cpp ../CTPlottingScripts/src/ctTrainingPlots.cpp)
add_executable(BuildDBParamCombinations ${SRC})
target_link_libraries(BuildDBParamCombinations patterns rootTreeFiles)

# Builds all template banks (template dbs) given inside the datapoints vector in configuration.h
set(SRC1 main_CTbuildDatapoints.cpp src/cosmicTemplatesBuild.cpp ../CTPlottingScripts/src/ctTrainingPlots.cpp)
add_executable(BuildDBMultiConfigs ${SRC1})
target_link_libraries(BuildDBMultiConfigs patterns rootTreeFiles)

# builds a single template bank (template db) according to command line parameters
set(SRC2 main_CTbuild.cpp src/cosmicTemplatesBuild.cpp ../root/TemplateDatabaseFile.cpp)
add_executable(BuildDBSingleConfig ${SRC2})
target_link_libraries(BuildDBSingleConfig patterns rootTreeFiles)

# evaluates the cosmic efficiency of an already build template db using a separate mc cosmic dataset
set(SRC3 main_CosmicEffTemplFilter.cpp src/cosmicEffTemplFilter.cpp src/getCosmicSIDtracks.cpp ../root/TemplateDatabaseFile.cpp)
add_executable(CosmicEffTemplFilter ${SRC3})
target_link_libraries(CosmicEffTemplFilter patterns rootTreeFiles)

# slims down mc cosmic datasets from slim segs (from mu3e simulation/trirec)
set(SRC5 main_getCosmicSIDtracks.cpp src/cosmicEffTemplFilter.cpp src/getCosmicSIDtracks.cpp)
add_executable(MakeCosmicSIDTrackFile ${SRC5})
target_link_libraries(MakeCosmicSIDTrackFile rootTreeFiles patterns)