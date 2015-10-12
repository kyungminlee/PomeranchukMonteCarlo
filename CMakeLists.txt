cmake_minimum_required(VERSION 3.2)
project(pomeranchuk)

#find_package(LAPACK REQUIRED)

include_directories(
        /Users/kmlee/.brew/include
        /Users/kmlee/.brew/include/eigen3
        /Users/kmlee/.brew/opt/openblas/include
        ../koretools/include
)

link_directories(
        /Users/kmlee/.brew/opt/openblas/lib
)

add_compile_options(-std=c++11)

link_libraries(openblas)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    src/main.cc
    src/pome/basics.h
    src/pome/geometry.h
    src/pome/geometry_impl.h
    src/pome/measurement.h
    src/pome/particle.h
    src/pome/pomeranchuk.h
    src/pome/theta.h
    src/pome/wavefunction.h
)


add_executable(pomeranchuk ${SOURCE_FILES})