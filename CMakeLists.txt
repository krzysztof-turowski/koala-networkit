cmake_minimum_required(VERSION 3.20)

find_program(CMAKE_CXX_COMPILER NAMES g++ clang++)

project(koala-networkit LANGUAGES CXX)
set(PROJECT_NAME koala-networkit)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(KOALA_CXX_FLAGS "-Wall -Wextra -Wstrict-aliasing -Wpedantic -Werror -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-overflow=2 -Wswitch-default -Wundef -Wno-error -fpermissive -fdiagnostics-show-option -O3 -pthread")
set(KOALA_LINKER_FLAGS "")
message("Compiled with ${CMAKE_CXX_COMPILER_ID}: ${CMAKE_CXX_COMPILER}")
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(KOALA_CXX_FLAGS "${KOALA_CXX_FLAGS} -ferror-limit=5")
    set(KOALA_LINKER_FLAGS "${KOALA_LINKER_FLAGS} -Xclang -fopenmp")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(KOALA_CXX_FLAGS "${KOALA_CXX_FLAGS} -fmax-errors=5 -Wlogical-op -Wnoexcept -Wno-unknown-pragmas -Wstrict-null-sentinel -lstdc++")
    set(KOALA_LINKER_FLAGS "${KOALA_LINKER_FLAGS} -fopenmp")
endif()

set(CMAKE_CXX_CPPLINT cpplint --linelength=100 --extensions=cpp,hpp --filter=-legal/copyright,-build/c++20,-build/include,-build/namespaces,-build/header_guard,-runtime/references,-runtime/string --recursive --exclude=${CMAKE_SOURCE_DIR}/backup --exclude=${CMAKE_SOURCE_DIR}/build --exclude=${CMAKE_SOURCE_DIR}/input --exclude=${CMAKE_SOURCE_DIR}/lib --exclude=${CMAKE_SOURCE_DIR}/cpp/flow/maximum_flow/dynamic_tree/*)

include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

#<gtest>
FetchContent_Declare(
    gtest
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    GIT_TAG main
    GIT_PROGRESS TRUE
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/lib/gtest"
)
FetchContent_MakeAvailable(gtest)
#</gtest>

#<networKit>
FetchContent_Declare(
    networkit
    GIT_REPOSITORY "https://github.com/networkit/networkit.git"
    GIT_TAG master
    GIT_PROGRESS TRUE
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/lib/networkit"
)
FetchContent_MakeAvailable(networkit)
#</networKit>

#<boost>
FetchContent_Declare(
    boost
    GIT_REPOSITORY "https://github.com/boostorg/boost.git"
    GIT_TAG master
    GIT_PROGRESS TRUE
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/lib/boost"
)
FetchContent_MakeAvailable(boost)
#</boost>

#<Csdp>
FetchContent_Declare(
    csdp
    GIT_REPOSITORY "https://github.com/coin-or/Csdp.git"
    GIT_TAG master
    GIT_PROGRESS TRUE
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/lib/csdp"
)
FetchContent_MakeAvailable(csdp)
#</Csdp>

set(CMAKE_VERBOSE_MAKEFILE OFF)

add_custom_target(csdp ALL
    COMMAND CFLAGS=-I"${CMAKE_SOURCE_DIR}/lib/csdp/include" make libsdp.a
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/lib/csdp/lib")
add_library(csdp_lib STATIC IMPORTED)
set_target_properties(csdp_lib PROPERTIES
    IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/csdp/lib/libsdp.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/lib/csdp/include"
    LINK_FLAGS ${KOALA_LINKER_FLAGS})
add_dependencies(csdp_lib csdp)

include_directories("${CMAKE_SOURCE_DIR}/include")

option(KOALA_BUILD_SHARED "Build shared library" OFF)
option(KOALA_BUILD_TESTS "Build tests" ON)
option(KOALA_BUILD_BENCHMARKS "Build benchmarks" ON)

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
target_link_libraries(${PROJECT_NAME} networkit boost_dynamic_bitset boost_graph csdp_lib lapack blas gfortran quadmath)
set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS ${KOALA_CXX_FLAGS}
    LINK_FLAGS ${KOALA_LINKER_FLAGS})
add_subdirectory(cpp)

function(koala_make_test testname testfile)
    add_executable(${testname} ${CMAKE_CURRENT_LIST_DIR}/${testfile})
    target_link_libraries(${testname} ${PROJECT_NAME} gtest_main)
    set_target_properties(${testname} PROPERTIES
        COMPILE_FLAGS ${KOALA_CXX_FLAGS}
        LINK_FLAGS ${KOALA_LINKER_FLAGS})
    add_test(${testname} ${testname})
endfunction()

function(koala_make_benchmark testname testfile scriptfile)
    add_executable(${testname} ${CMAKE_CURRENT_LIST_DIR}/${testfile})
    target_link_libraries(${testname} ${PROJECT_NAME})
    set_target_properties(${testname} PROPERTIES
        COMPILE_FLAGS ${KOALA_CXX_FLAGS}
        LINK_FLAGS ${KOALA_LINKER_FLAGS})
    add_custom_command(TARGET ${testname}
        POST_BUILD COMMAND /bin/sh ${CMAKE_CURRENT_LIST_DIR}/${scriptfile} ${ARGN})
endfunction()

if(KOALA_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

if(KOALA_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()
