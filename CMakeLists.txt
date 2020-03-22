cmake_minimum_required(VERSION 3.15)

SET(PROJECT_ROOT_PATH "${CMAKE_SOURCE_DIR}/../")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/")
INCLUDE_DIRECTORIES("${PROJECT_ROOT_PATH}/include/")

project(NumericalAnalysis)

set(CMAKE_CXX_STANDARD 17)

add_executable(NumericalAnalysis main.cpp include/Chapter2/LinearEquations.h include/Type/Type.h)