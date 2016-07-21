cmake_minimum_required(VERSION 3.5)
project(test_1)

find_package(Boost 1.60.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
include_directories (include)
link_directories (lib)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SOURCE_FILES main.cpp)
add_executable(test_1 ${SOURCE_FILES})
target_link_libraries(test_1 ${Boost_LIBRARIES} bamtools)

