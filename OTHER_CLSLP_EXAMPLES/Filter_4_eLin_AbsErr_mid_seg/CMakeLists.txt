cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
# set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Set some basic project attributes
project ( CLSLP
          VERSION  1.00
          DESCRIPTION "CLSLP Linear-Phase FIR Filter Synthesis" )


file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)

# This project will output an executable file
add_executable(${PROJECT_NAME} ${SRC_FILES})

# Include the configuration header in the build
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")

# Math includes
target_link_libraries(${PROJECT_NAME} PUBLIC "gsl")
target_link_libraries(${PROJECT_NAME} PUBLIC "gslcblas")
target_link_libraries(${PROJECT_NAME} PUBLIC "fftw3")
target_link_libraries(${PROJECT_NAME} PUBLIC "fftw3_threads")

