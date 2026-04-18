cmake_minimum_required(VERSION 3.18)
project(leafblower_python CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 2.11 REQUIRED)

# Sources: all src/*.cpp EXCEPT r_bridge.cpp (includes Rinternals.h unavailable here)
set(CORE_SOURCES
    ../src/c_api.cpp
    ../src/logit.cpp
    ../src/lbfgsb_solver.cpp
    ../src/ieppa.cpp
)

pybind11_add_module(_leafblower
    leafblower/_bindings.cpp
    ${CORE_SOURCES}
)

target_include_directories(_leafblower PRIVATE ../src)
target_compile_options(_leafblower PRIVATE -O3)
