cmake_minimum_required(VERSION 3.20)
project(koala-networkit CXX)
find_program(MAKE_EXE NAMES make gmake nmake)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(LIB_DIR "${CMAKE_SOURCE_DIR}/lib")

include(FetchContent)

#<gtest>
FetchContent_Declare(
        gtest
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        GIT_TAG master
        SOURCE_DIR "${LIB_DIR}/gtest"
)
FetchContent_MakeAvailable(gtest)
#</gtest>

#<networKit>
FetchContent_Declare(
        networkit
        GIT_REPOSITORY "https://github.com/networkit/networkit.git"
        GIT_TAG master
        SOURCE_DIR "${LIB_DIR}/networkit"
)
FetchContent_MakeAvailable(networkit)
#</networKit>

option(KOALA_BUILD_CORE "Build Koala library" ON)
option(KOALA_BUILD_TESTS "Build Koala tests" ON)

if (KOALA_BUILD_CORE)
    add_library(koala-networkit STATIC cpp/koala.cpp)
    target_include_directories(koala-networkit PUBLIC cpp/io)
    target_include_directories(koala-networkit PUBLIC include)
    target_link_libraries(koala-networkit networkit)
endif()

if (KOALA_BUILD_TESTS)
    enable_testing()

    add_executable(test_graph_io test/testGraphIO.cpp)
    target_include_directories(test_graph_io PUBLIC include)
    target_link_libraries(test_graph_io gtest_main networkit)
    add_test(test_graph_io test_graph_io)
endif()
