cmake_minimum_required(VERSION 3.20)
project(koala-networkit CXX)
set(PROJECT_NAME "koala-networkit")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(FetchContent)

#<gtest>
FetchContent_Declare(
        gtest
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        GIT_TAG master
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/lib/gtest"
)
FetchContent_MakeAvailable(gtest)
#</gtest>

#<networKit>
FetchContent_Declare(
        networkit
        GIT_REPOSITORY "https://github.com/networkit/networkit.git"
        GIT_TAG master
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/lib/networkit"
)
FetchContent_MakeAvailable(networkit)
#</networKit>

option(KOALA_BUILD_SHARED "Build shared library" OFF)
option(KOALA_BUILD_TESTS "Build tests" ON)

function(koala_add_module modname)
    foreach(file ${ARGN})
        target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/${file})
    endforeach()
endfunction()

if(KOALA_BUILD_SHARED)
    add_library(${PROJECT_NAME} SHARED cpp/koala.cpp)
else()
    add_library(${PROJECT_NAME} STATIC cpp/koala.cpp)
endif()
target_include_directories(${PROJECT_NAME} PUBLIC include)
target_link_libraries(${PROJECT_NAME} networkit)
add_subdirectory(cpp)

function(koala_make_test testname testfile)
    add_executable(${testname} ${CMAKE_CURRENT_LIST_DIR}/${testfile})
    target_include_directories(${testname} PUBLIC include)
    target_link_libraries(${testname} gtest_main networkit ${PROJECT_NAME})
    add_test(${testname} ${testname})
endfunction()

if(KOALA_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
