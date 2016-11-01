cmake_minimum_required(VERSION 3.5)
project(geep)

find_package(Boost 1.60.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
include_directories (include)
link_directories (lib)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SOURCE_FILES src/main.cpp src/bam_reader.cpp include/bam_reader.h src/annotation_reader.cpp include/annotation_reader.h src/interval_map.cpp include/interval_map.h include/string_tools.h include/rpkm_calculation.h src/rpkm_calculation.cpp src/test.cpp include/test.h)

add_executable(geep ${SOURCE_FILES})
target_link_libraries(geep ${Boost_LIBRARIES} bamtools)

