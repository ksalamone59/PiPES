cmake_minimum_required(VERSION 3.10)
project(PiPES)
add_library(PiPES INTERFACE)
target_include_directories(PiPES INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_definitions(PiPES INTERFACE
    PHASE_SHIFT_DATADIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
)
target_compile_features(PiPES INTERFACE cxx_std_20)