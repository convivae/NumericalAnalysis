cmake_minimum_required(VERSION 3.15)

SET(PROJECT_ROOT_PATH "${CMAKE_SOURCE_DIR}/../")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/")

project(NumericalAnalysis)

set(CMAKE_CXX_STANDARD 17)

add_executable(NumericalAnalysis main.cpp
        Chapter2/LinearEquations.h
        Type/Type.h
        Chapter3/MatrixEigenvaluesAndEigenvectors.cpp
        Chapter3/MatrixEigenvaluesAndEigenvectors.h
        lib/MatrixOperation.h)