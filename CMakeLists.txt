cmake_minimum_required(VERSION 3.20)

find_program(CMAKE_CXX_COMPILER NAMES g++ clang++)

project(koala-networkit LANGUAGES CXX)
set(PROJECT_NAME koala-networkit)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(KOALA_CXX_FLAGS "-Wall -Wextra -Wstrict-aliasing -Wpedantic -Wno-error -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-overflow=2 -Wswitch-default -Wundef -fpermissive -fdiagnostics-show-option -pthread")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(KOALA_CXX_FLAGS "${KOALA_CXX_FLAGS} -g -O0")
else()
    set(KOALA_CXX_FLAGS "${KOALA_CXX_FLAGS} -O3")
endif()
set(KOALA_LINKER_FLAGS "")

message("Compiled with ${CMAKE_CXX_COMPILER_ID}: ${CMAKE_CXX_COMPILER}")
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(KOALA_CXX_FLAGS "${KOALA_CXX_FLAGS} -ferror-limit=5")
    set(KOALA_LINKER_FLAGS "${KOALA_LINKER_FLAGS} -Xclang -fopenmp")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(KOALA_CXX_FLAGS "${KOALA_CXX_FLAGS} -fmax-errors=5 -Wlogical-op -Wnoexcept -Wno-unknown-pragmas -Wstrict-null-sentinel -lstdc++")
    set(KOALA_LINKER_FLAGS "${KOALA_LINKER_FLAGS} -fopenmp")
endif()

set(CMAKE_CXX_CPPLINT cpplint --quiet --linelength=100 --extensions=cpp,hpp --filter=-legal/copyright,-build/c++11,-build/c++17,-build/c++20,-runtime/references,-runtime/string --recursive --exclude=${CMAKE_SOURCE_DIR}/backup --exclude=${CMAKE_SOURCE_DIR}/build --exclude=${CMAKE_SOURCE_DIR}/input --exclude=${CMAKE_SOURCE_DIR}/lib)

include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

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

set(NETWORKIT_SOURCE_DIRS
  "${CMAKE_SOURCE_DIR}/lib/networkit/include"
  "${CMAKE_SOURCE_DIR}/lib/networkit/extlibs/dhb/include"
  "${CMAKE_SOURCE_DIR}/lib/networkit/extlibs/ttmath"
  "${CMAKE_SOURCE_DIR}/lib/networkit/extlibs/tlx"
)

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
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_custom_target(csdp ALL
    COMMAND CFLAGS=-I"${CMAKE_SOURCE_DIR}/lib/csdp/include" make libsdp.a
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/lib/csdp/lib")
add_library(csdp_lib STATIC IMPORTED)
set_target_properties(csdp_lib PROPERTIES
    IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/csdp/lib/libsdp.a"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/lib/csdp/include"
    LINK_FLAGS ${KOALA_LINKER_FLAGS})
add_dependencies(csdp_lib csdp)

option(KOALA_BUILD_SHARED "Build shared library" OFF)
option(KOALA_BUILD_TESTS "Build tests" ON)
option(KOALA_BUILD_BENCHMARKS "Build benchmarks" ON)

function(koala_add_module modname)
    foreach(file ${ARGN})
        target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${file})
    endforeach()
endfunction()

function(koala_common_options target)
    target_include_directories(${target} PUBLIC "${CMAKE_SOURCE_DIR}/include")
    target_include_directories(${target} SYSTEM PRIVATE ${NETWORKIT_SOURCE_DIRS})
    set_target_properties(${target} PROPERTIES
        COMPILE_FLAGS ${KOALA_CXX_FLAGS}
        LINK_FLAGS ${KOALA_LINKER_FLAGS})
endfunction()

if(KOALA_BUILD_SHARED)
    add_library(${PROJECT_NAME} SHARED cpp/koala.cpp)
else()
    add_library(${PROJECT_NAME} STATIC cpp/koala.cpp)
endif()
koala_common_options(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}
    networkit ntl Eigen3::Eigen boost_graph csdp_lib lapack blas gfortran quadmath)
add_subdirectory(cpp)

function(koala_make_test testname testfile)
    add_test(NAME ${testname} COMMAND ${testname})
    set_tests_properties(${testname} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
    add_executable(${testname} ${testfile})
    koala_common_options(${testname})
    target_include_directories(${testname} PRIVATE "${CMAKE_SOURCE_DIR}")
    target_link_libraries(${testname} ${PROJECT_NAME} gtest gtest_main)
endfunction()

function(koala_make_benchmark testname testfile scriptfile)
    add_executable(${testname} ${testfile})
    koala_common_options(${testname})
    target_link_libraries(${testname} ${PROJECT_NAME})
    add_custom_command(TARGET ${testname}
        POST_BUILD COMMAND /bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/${scriptfile} ${ARGN})
endfunction()

if(KOALA_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

if(KOALA_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()
