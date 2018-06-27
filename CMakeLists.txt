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

set(MECACELL_DIR include/mecacell)
set(MECACELL_VIEWER_DIR ${MECACELL_DIR}/mecacellviewer)

add_executable(console ${SRC} src/mainconsole.cpp)
target_include_directories(console PUBLIC ${MECACELL_DIR})

# add_executable(spherepacking ${SRC} src/mainspherepacking.cpp)
# target_link_libraries(spherepacking mecacell)

find_package(Qt5 COMPONENTS Core Quick REQUIRED)
set(CMAKE_AUTOMOC ON)
cmake_policy(SET CMP0071 NEW)
qt5_wrap_cpp(MOC_SRC ${MECACELL_VIEWER_DIR}/signalslotbase.h ${MECACELL_VIEWER_DIR}/camera.hpp)
qt5_add_resources(RESOURCES ${MECACELL_VIEWER_DIR}/resourcesLibMecacellViewer.qrc)

add_executable(viewer ${SRC} ${MOC_SRC} ${RESOURCES} src/mainviewer.cpp)
target_link_libraries(viewer Qt5::Core Qt5::Quick)
target_include_directories(viewer PUBLIC ${MECACELL_DIR} ${MECACELL_VIEWER_DIR})
