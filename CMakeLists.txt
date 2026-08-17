cmake_minimum_required(VERSION 3.12)

project(PyCARDModules)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -fPIC")

# -------------------------------------------------------------------
# Dependencies
# -------------------------------------------------------------------

find_package(Python COMPONENTS Development Interpreter REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Armadillo REQUIRED)

# -------------------------------------------------------------------
# CARD reference-based module
# -------------------------------------------------------------------

add_library(CARDref_module MODULE
    src/CARDref.cpp
    src/CARDref_wrapper.cpp
)

target_include_directories(CARDref_module PRIVATE
    ${pybind11_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(CARDref_module PRIVATE
    ${ARMADILLO_LIBRARIES}
)

set_target_properties(CARDref_module PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)

# -------------------------------------------------------------------
# CARD reference-free module
# -------------------------------------------------------------------

add_library(CARDfree_module MODULE
    src/CARDfree.cpp
    src/CARDfree_wrapper.cpp
)

target_include_directories(CARDfree_module PRIVATE
    ${pybind11_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(CARDfree_module PRIVATE
    ${ARMADILLO_LIBRARIES}
)

set_target_properties(CARDfree_module PROPERTIES
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)

# -------------------------------------------------------------------
# Platform-specific settings
# -------------------------------------------------------------------

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

if(UNIX AND NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -shared")

    target_link_libraries(CARDref_module PRIVATE
        ${CMAKE_DL_LIBS}
    )

    target_link_libraries(CARDfree_module PRIVATE
        ${CMAKE_DL_LIBS}
    )
endif()