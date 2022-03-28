cmake_minimum_required(VERSION 3.21)
project(lab2)

set(CMAKE_CXX_STANDARD 23)

include_directories(../software)

add_executable(lab1 main.cpp ../software/DECOMP.CPP)