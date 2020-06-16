CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

add_executable(Mu3eCosPat PatternEngine.cpp main.cpp SuperPixel.cpp SuperPixel.h)