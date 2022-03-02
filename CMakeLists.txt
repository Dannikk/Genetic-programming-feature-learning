cmake_minimum_required(VERSION 3.20)
project(GPFL)

#set(CMAKE_CXX_STANDARD 14)

add_definitions(-DNO_DLL)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99 -fopenmp")

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(GPFL main.cpp cgp.c)

target_link_libraries( GPFL ${OpenCV_LIBS} )

#message(${OpenCV_INCLUDE_DIRS})
#message(${OpenCV_LIBS})

#message(STATUS "Checking OpenMP")
#find_package(OpenMP)
#IF(OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
##    message("Found OpenMP!")
#endif()