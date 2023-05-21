project(banana)
cmake_minimum_required(VERSION 3.8)

# put binary output in ./ directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(LIBRARY_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# select release build type as default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "'${CMAKE_BUILD_TYPE}' built type was selected.")

# add include directories
include_directories(SYSTEM "${PROJECT_SOURCE_DIR}")

# add compiler flags to define C version
# and to compile position independent code
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c89 -fPIC")

# add compiler flags to enable strict compilation in debug mode
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 --coverage")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3 ")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra -pedantic")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wmissing-prototypes")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wstrict-prototypes")

# add compiler flags to define C++ version
# and to compile position independent code
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

# add C++ compiler flag to enable code coverage testing
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")

set(HEADERS
    ${PROJECT_SOURCE_DIR}/src/util/lodepng.h
    ${PROJECT_SOURCE_DIR}/src/catch.hpp
    ${PROJECT_SOURCE_DIR}/src/util/fitsfile.hpp
    ${PROJECT_SOURCE_DIR}/src/hedgehog.hpp
    ${PROJECT_SOURCE_DIR}/src/originalBanana.hpp
    ${PROJECT_SOURCE_DIR}/src/papaya2.hpp
    ${PROJECT_SOURCE_DIR}/src/peach.hpp
    ${PROJECT_SOURCE_DIR}/src/photofunctions.hpp
    ${PROJECT_SOURCE_DIR}/src/readarg.hpp
    ${PROJECT_SOURCE_DIR}/src/readdata.hpp
    ${PROJECT_SOURCE_DIR}/src/util/savepng.hpp
    ${PROJECT_SOURCE_DIR}/src/sugarsnappea.hpp
    ${PROJECT_SOURCE_DIR}/src/tools.hpp
    ${PROJECT_SOURCE_DIR}/src/vector.hpp
)

# add source files of this module
set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/util/savepng.cpp
    ${PROJECT_SOURCE_DIR}/src/util/lodepng.cpp
    ${PROJECT_SOURCE_DIR}/src/util/fitsfile.cpp
    ${PROJECT_SOURCE_DIR}/src/banana.cpp
)

# add object files to be compiled
add_executable(banana ${SOURCES})
#target_link_libraries(banana "stdc++fs")# "c++experimental")
target_link_libraries(banana cfitsio)
target_link_libraries(banana CCfits)

