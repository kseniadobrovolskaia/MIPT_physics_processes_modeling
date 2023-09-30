cmake_minimum_required(VERSION 3.14)

project(HarmonicSimulator)

add_compile_options(-Wall -std=c++2a -fexceptions)

set(SOURCE_EXE ./HarmonicSimulator/main.cpp)			

include_directories(include)			

add_executable(Simulator ${SOURCE_EXE})

add_subdirectory(definitions)				
#add_subdirectory(Tests)

target_link_libraries(Simulator HarmonicSimulator)
