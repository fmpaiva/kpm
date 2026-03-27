cmake_minimum_required(VERSION 3.31)
project(kpm)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

include_directories(.)

add_executable(kpm
        hamiltonian.h
        kpm.cpp
        kpm.h
        main.cpp
        random.h
        simulation.cpp
        simulation.h
        Constants.h
        hamiltonian.cpp
)
