cmake_minimum_required(VERSION 3.10)
project(ModuleABM)

set(CMAKE_CXX_STANDARD 14)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

set(BOOST_ROOT C:/local/boost_1_64_0/)
set(BOOST_LIBRARYDIR ${BOOST_ROOT}/lib64-msvc-14.1/)

find_package(Boost 1.64 REQUIRED COMPONENTS system thread)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(ModuleABM main.cpp structs.h functions.h)
    target_link_libraries(ModuleABM ${Boost_LIBRARIES})
endif()
