CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)



add_library(patterns SHARED SuperPixelHit.cpp PatternEngine.cpp)
target_link_libraries(patterns ${ROOT_LIBRARIES})
#link_directories(${ROOT_LIBRARY_DIR})

add_executable(Mu3eCosPat main.cpp)
target_link_libraries(Mu3eCosPat patterns)