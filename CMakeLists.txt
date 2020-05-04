cmake_minimum_required(VERSION 3.15)
project(MPI)

set(CMAKE_CXX_STANDARD 17)



find_package(MPI REQUIRED)
# add this line only when you are using openmpi which has a different c++ bindings
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Use imported targets would make things much eazier. Thanks Levi for pointing it out.
add_executable(MPI main.cpp)

target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES})