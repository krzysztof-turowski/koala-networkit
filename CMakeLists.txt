cmake_minimum_required(VERSION 3.20)
project(koala-networkit LANGUAGES CXX)
set(PROJECT_NAME koala-networkit)

find_program(CMAKE_CXX_COMPILER NAMES g++ clang++)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(KOALA_CXX_FLAGS "-Wall -Wextra -Wstrict-aliasing -Wpedantic -Werror -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-overflow=2 -Wswitch-default -Wundef -Wno-error -fpermissive -fdiagnostics-show-option -O3 -pthread")
set(KOALA_LINKER_FLAGS "")
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(KOALA_CXX_FLAGS "${KOALA_CXX_FLAGS} -ferror-limit=5")
    set(KOALA_LINKER_FLAGS "${KOALA_LINKER_FLAGS} -Xclang -fopenmp -lomp")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(KOALA_CXX_FLAGS "${KOALA_CXX_FLAGS} -fmax-errors=5 -Wlogical-op -Wnoexcept -Wstrict-null-sentinel -lstdc++")
    set(KOALA_LINKER_FLAGS "${KOALA_LINKER_FLAGS} -fopenmp")
endif()

set(CMAKE_CXX_CPPLINT cpplint --linelength=100 --extensions=cpp,hpp --filter=-legal/copyright,-build/c++20,-build/include,-build/namespaces,-build/header_guard,-runtime/references,-runtime/string --recursive)
set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,readability-*)

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
set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS ${KOALA_CXX_FLAGS}
    LINK_FLAGS ${KOALA_LINKER_FLAGS})
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
