cmake_minimum_required(VERSION 2.8)
project(lala)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "-O3 -std=c++11 -Wall -Wextra -pedantic")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

file(GLOB SRC
	"src/*.h"
	"src/*.hpp"
)

add_executable(lala ${SRC} src/main.cpp)
