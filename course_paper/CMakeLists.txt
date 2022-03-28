cmake_minimum_required(VERSION 3.21)
project(coursework)

set(CMAKE_CXX_STANDARD 23)

include_directories(../software)

add_executable(coursework main.cpp ../software/ZEROIN.CPP ../software/QUANC8.CPP ../software/RKF45.CPP)
