cmake_policy(SET CMP0048 NEW)
cmake_minimum_required(VERSION 3.0)
#set(CMAKE_CXX_STANDARD 11)
set (CMAKE_C_COMPILER "/usr/bin/gcc")
set (CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_FLAGS -std=c++11)
#add_compile_options(-std=c++11 -Wall -Werror)

#get_filename_component(example_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(xx)
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build)

# Let's require mutation++. It should export a suitable config that also
# exports include directories for Eigen. The `if` is needed to allow the
# examples to build during the process of building the entire mutation++
# tree.
if(NOT TARGET mutation++)
    find_package(mutation++ REQUIRED)
endif()
find_package(Eigen3 REQUIRED)
#message("Eigen++ = ${Eigen3_FOUND},include= ${EIGEN3_INCLUDE_DIRS}")
#message("mutation++ = ${mutation++_FOUND},include= ${mutation++_LIBRARIES}")
include_directories($ENV{MPP_DIRECTORY}/install/include/mutation++)
aux_source_directory(src/main DIR_MAIN)
aux_source_directory(src/sub DIR_SUB)
add_executable(xx ${DIR_MAIN})
include_directories("${PROJECT_SOURCE_DIR}/src/sub")
add_library (sub STATIC ${DIR_SUB})
target_link_libraries(xx
    PRIVATE
        mutation++
        sub
)
