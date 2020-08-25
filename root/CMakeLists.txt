CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(mu3eRecAcc C CXX)
set(CMAKE_CXX_STANDARD 14)

add_library(rootTreeFiles STATIC SegsTreeRead.cpp SlimSegsTree.cpp Mu3eTree.cpp Mu3eTree.h)
target_link_libraries(rootTreeFiles ${ROOT_LIBRARIES})