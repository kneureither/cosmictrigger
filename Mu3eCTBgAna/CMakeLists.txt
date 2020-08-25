CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(inc)
include_directories(src)

set(SRC main.cpp)
add_executable(Mu3eCTBgAna ${SRC})
target_link_libraries(Mu3eCTBgAna patterns rootTreeFiles)