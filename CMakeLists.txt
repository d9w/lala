cmake_minimum_required(VERSION 2.8)
project(lala)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "-O3 -pthread -std=c++14 -pedantic")
set(CMAKE_CXX_STANDARD 14)
set(MecaCellDir ~/MecaCell)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

file(GLOB SRC
	"src/scenario.hpp"
	"src/cell.hpp"
	"src/config.hpp"
	)

link_directories(${MecaCellDir}/build/mecacell)
include_directories(${MecaCellDir})

add_executable(console ${SRC} src/mainconsole.cpp)
target_link_libraries(console mecacell)

add_executable(spherepacking ${SRC} src/mainspherepacking.cpp)
target_link_libraries(spherepacking mecacell)

set(CMAKE_AUTOMOC ON)
link_directories(${MecaCellDir}/build/mecacellviewer)
find_package(Qt5Core)
find_package(Qt5Gui)
find_package(Qt5Quick)
find_package(Qt5OpenGL)
add_executable(viewer ${SRC} src/mainviewer.cpp)
target_link_libraries(viewer mecacell)
target_link_libraries(viewer mecacellviewer)
qt5_use_modules(viewer Quick Core Gui OpenGL)
