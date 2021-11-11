cmake_minimum_required(VERSION 3.5)
project(bamsplice)

set( EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" )
find_package(Boost COMPONENTS system filesystem thread REQUIRED)

add_subdirectory(./3rd_party/bamtools)
include_directories(${Boost_INCLUDE_DIRS})
include_directories (include)
link_directories (lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SOURCE_FILES src/main.cpp
                 src/bam_reader.cpp
                 include/bam_reader.h
                 src/annotation_reader.cpp
                 include/annotation_reader.h
                 src/interval_map.cpp
                 include/interval_map.h
                 include/string_tools.h
                 include/rpkm_calculation.h
                 src/rpkm_calculation.cpp
                 src/test.cpp
                 include/test.h
                 src/thread.cpp
                 include/thread.h
                 include/cxxopts/cxxopts.hpp)

add_executable(bamsplice ${SOURCE_FILES})
add_dependencies (bamsplice BamTools)
target_link_libraries(bamsplice ${Boost_LIBRARIES} bamtools)

install(TARGETS bamsplice
        RUNTIME DESTINATION bin)