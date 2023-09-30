cmake_minimum_required(VERSION 3.14)

project(definitions)

add_compile_options(-Wall -std=c++2a -fexceptions)

include_directories(.././include)			

set(SOURCE_LIB  Solver.cpp
				DiffEquation.cpp

	)

add_library(HarmonicSimulator STATIC ${SOURCE_LIB})


