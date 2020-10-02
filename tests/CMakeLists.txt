CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(mu3eRecAcc C CXX)
set(CMAKE_CXX_STANDARD 14)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/root)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Mu3eCosPat/)

set(SRC tests.cpp ../root/SlimSegsTree.cpp testTemplateID.h ../Mu3eCosPat/TemplateData.cpp)
add_executable(tests ${SRC})
target_link_libraries(tests ${ROOT_LIBRARIES} f2clib karilib patterns)
add_compile_options("-lf2c" "-DMAIN")
#target_compile_options(Mu3eRecAcc PUBLIC -lf2c -DMAIN)

add_executable(testRootMCtype findRootMctypeError.cpp ../root/SegsTreeRead.cpp)
target_link_libraries(testRootMCtype ${ROOT_LIBRARIES})
