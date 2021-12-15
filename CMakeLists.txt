cmake_minimum_required(VERSION 3.20)
project(LDScvx)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS "-O3")

add_executable(LDScvx main.cpp Graph.cpp Graph.h FlowNetwork.cpp FlowNetwork.h EdgeFN.cpp EdgeFN.h Args.cpp Args.h)
