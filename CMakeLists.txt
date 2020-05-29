cmake_minimum_required(VERSION 3.15)

SET(PROJECT_ROOT_PATH "${CMAKE_SOURCE_DIR}/../")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/")

project(NumericalAnalysis)

set(CMAKE_CXX_STANDARD 17)

add_executable(NumericalAnalysis
        main.cpp
        Type/Type.h
        lib/MatrixOperation.h
        lib/fun.h
        lib/Polymerization.h
        lib/Polymerization.cpp
        Chapter2/LinearEquations.h
        Chapter3/MatrixEigenvaluesAndEigenvectors.cpp
        Chapter3/MatrixEigenvaluesAndEigenvectors.h
        Chapter4/NonlinearEquations.cpp
        Chapter4/NonlinearEquations.h
        Chapter5/InterpolationAndApproximation.cpp
        Chapter5/InterpolationAndApproximation.h
        Chapter6/NumericalIntegration.cpp
        Chapter6/NumericalIntegration.h
        lib/OrthogonalPolynomial.cpp
        lib/OrthogonalPolynomial.h
        )