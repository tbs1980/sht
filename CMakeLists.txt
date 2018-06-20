cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(sht VERSION 0.1.0 LANGUAGES CXX)

# vscode cmake helper
include(CMakeToolsHelpers OPTIONAL)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(include)
add_subdirectory(test)
add_subdirectory(bench)
