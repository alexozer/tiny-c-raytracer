cmake_minimum_required(VERSION 3.12)
project(raytrace C)

set(CMAKE_C_STANDARD 11)

add_executable(raytrace src/main.c src/export.c src/export.h)
target_link_libraries(raytrace m png)