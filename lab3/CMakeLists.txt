cmake_minimum_required(VERSION 3.21)
project(lab3)

set(CMAKE_CXX_STANDARD 23)

include_directories(../software)

add_executable(lab3 main.cpp ../software/RKF45.CPP ../software/FORSYTHE.H)
