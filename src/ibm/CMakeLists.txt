# this file is used for the programme CMake,
# which allows us to build and compile the executable
# see https://cmake.org/cmake/help/latest/index.html
cmake_minimum_required(VERSION 3.10)


#### VARIABLES ####

# put executable name in variable, as we
# need it in multiple places
set(EXE_NAME_IND_TEST "sexsel_space_ind_tests")
set(EXE_NAME_SIM_TEST "sexsel_space_sim_tests")
set(EXE_NAME_PATCH_TEST "sexsel_space_patch_tests")
set(EXE_NAME "sexsel_space")


# google test url
set(GOOGLETEST_URL "https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip") 

#### THE ACTUAL PROJECT FILE ####

# set the project name and version
project(sexsel_space_multisignal)

add_library(sexsel_library individual.cpp simulation.cpp patch.cpp)

# make sure it is c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# now adding the google mock testing framework
# see: https://google.github.io/googletest/quickstart-cmake.html

# call cmake's FetchContent module
# see https://cmake.org/cmake/help/latest/module/FetchContent.html 
include(FetchContent)
FetchContent_Declare(
    googletest
    URL "${GOOGLETEST_URL}"
)

FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(
    "${EXE_NAME}"
    main.cpp
)

target_link_libraries(
    "${EXE_NAME}" 
    sexsel_library
)

target_include_directories(
    "${EXE_NAME}" 
    PUBLIC
    "${PROJECT_BINARY_DIR}"
    )


add_executable(
    "${EXE_NAME_IND_TEST}"
    individual_tests.cpp
)

target_link_libraries(
    "${EXE_NAME_IND_TEST}" 
    gtest_main
    sexsel_library
)

target_include_directories(
    "${EXE_NAME_IND_TEST}" 
    PUBLIC
    "${PROJECT_BINARY_DIR}"
    )

add_executable(
    "${EXE_NAME_SIM_TEST}"
    simulation_tests.cpp
)

target_link_libraries(
    "${EXE_NAME_SIM_TEST}" 
    gtest_main 
    sexsel_library
)


target_include_directories(
    "${EXE_NAME_SIM_TEST}" 
    PUBLIC
    "${PROJECT_BINARY_DIR}"
    )

include(GoogleTest)
gtest_discover_tests("${EXE_NAME_IND_TEST}")
gtest_discover_tests("${EXE_NAME_SIM_TEST}")

