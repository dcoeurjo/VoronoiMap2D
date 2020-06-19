cmake_minimum_required(VERSION 3.1)

project(VornoiMap2D VERSION 1.0
        DESCRIPTION "Fast Voronoi Map in 2D"
        LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
                
add_executable(example example)
