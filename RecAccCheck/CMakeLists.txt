CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(CosmicTrigger C CXX)
set(CMAKE_CXX_STANDARD 14)

set(SRC recacc_main.cpp reconstructionAccuracyScript.cpp reconstruction_accuracy.cpp ../root/SegsTreeRead.cpp)
add_executable(RecAccCheck ${SRC})
target_link_libraries(RecAccCheck ${ROOT_LIBRARIES} f2clib karilib)
add_compile_options("-lf2c" "-DMAIN")
#target_compile_options(RecAccCheck PUBLIC -lf2c -DMAIN)

