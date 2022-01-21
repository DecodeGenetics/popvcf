cmake_minimum_required (VERSION 3.2)
project (popvcf LANGUAGES CXX)

include(ExternalProject)

# Build popvcf in release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Choose the type of build." FORCE)
endif()

# Build paw (dependency) in release by default
if(NOT PAW_CMAKE_BUILD_TYPE)
  set(PAW_CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Choose the type of build for paw." FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

############
## popvcf ##
############
set (popvcf_VERSION_MAJOR 0)
set (popvcf_VERSION_MINOR 0)
set (popvcf_VERSION_PATCH 1)

add_subdirectory(src) # Exposes "popvcf_sources", which contains all source files of popvcf
add_library(popvcf_objects OBJECT ${popvcf_sources})

target_compile_features(popvcf_objects PUBLIC cxx_std_17)
target_compile_options(popvcf_objects PUBLIC -Wall -Wextra -Wfatal-errors -pedantic -Wno-variadic-macros)
add_dependencies(popvcf_objects project_paw) # popvcf depends on building paw

target_include_directories(popvcf_objects PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_BINARY_DIR}/include/popvcf)

# Add popvcf executable
add_executable(popvcf src/main.cpp $<TARGET_OBJECTS:popvcf_objects>)

# Inherit stuff from popvcf_objects
target_include_directories(popvcf PRIVATE $<TARGET_PROPERTY:popvcf_objects,INTERFACE_INCLUDE_DIRECTORIES>)
target_compile_features(popvcf PRIVATE $<TARGET_PROPERTY:popvcf_objects,COMPILE_FEATURES>)
target_compile_options(popvcf PRIVATE $<TARGET_PROPERTY:popvcf_objects,COMPILE_OPTIONS>)

# configure a header file to pass some of the CMake settings to the source code
configure_file (
  ${PROJECT_SOURCE_DIR}/src/in.constants.hpp
  ${PROJECT_BINARY_DIR}/include/popvcf/constants.hpp)

target_include_directories(popvcf_objects PUBLIC ${PROJECT_BINARY_DIR}/include)

########################
## Build dependencies ##
########################
### paw ###
if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  message (STATUS "paw in debug mode")

  ExternalProject_Add(
  project_paw
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw
  BUILD_IN_SOURCE 0
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/paw
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -H${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw -B${CMAKE_CURRENT_BINARY_DIR}/paw -DCMAKE_BUILD_TYPE=${PAW_CMAKE_BUILD_TYPE} -DFORCE_AVX_FAST=1
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/paw static
  INSTALL_COMMAND "")
else()
  message (STATUS "paw in non-debug mode")

  ExternalProject_Add(
  project_paw
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw
  BUILD_IN_SOURCE 0
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/paw
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -H${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw -B${CMAKE_CURRENT_BINARY_DIR}/paw -DCMAKE_BUILD_TYPE=${PAW_CMAKE_BUILD_TYPE} -DFORCE_AVX512=1
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/paw static
  INSTALL_COMMAND "")
endif()

add_library(paw STATIC IMPORTED)
target_include_directories(popvcf_objects SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules/paw/include)
target_include_directories(popvcf_objects SYSTEM PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/paw/include)
set(paw_location ${CMAKE_CURRENT_BINARY_DIR}/paw/lib/libpaw.a)
message(STATUS "paw target location is ${paw_location}")
set_property(TARGET paw PROPERTY IMPORTED_LOCATION ${paw_location})
add_dependencies(paw project_paw)
add_dependencies(popvcf_objects paw)
target_link_libraries(popvcf PUBLIC paw)

### parallel_hashmap ###
target_include_directories(popvcf_objects SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/submodules/parallel-hashmap)

### threads ###
find_package(Threads)
target_link_libraries(popvcf PUBLIC ${CMAKE_THREAD_LIBS_INIT})

### rt and filesystem ###
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Using GCC")
    target_link_libraries(popvcf PUBLIC "rt")
    target_link_libraries(popvcf PUBLIC "stdc++fs")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message(STATUS "Using Clang")
  else()
    message(WARNING "Unsupported compiler")
  endif ()

#############
## Testing ##
#############

# TODO

################
# clang-format #
################
# Add clang-format check target
find_program(CLANG_FORMAT "clang-format-12")

if (CLANG_FORMAT STREQUAL "" OR CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
  find_program(CLANG_FORMAT "clang-format")
endif()

if (NOT CLANG_FORMAT STREQUAL "" AND NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
  add_custom_target (check_format "find" "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/test" "-name" "\"*.hpp\"" "-o" "-name" "\"*.cpp\"" "-exec" "${CLANG_FORMAT}" "--style=file" "--dry-run" "--Werror" "{}" "+" COMMENT "Checking all files format with clang-format.")
endif()

#################
# Documentation #
#################
## Adds a target to generate API documentation with Doxygen
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()

  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  configure_file(${doxyfile_in} ${doxyfile} @ONLY)

  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)

  #install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION share/doc)
endif()

###########
## Other ##
###########
# Get the current working git branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest abbreviated SHA commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_SHORT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest SHA commit hash of the working branch
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_LONG_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest SHA commit hash of the working branch
execute_process(
  COMMAND git diff COMMAND wc -l
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_NUM_DIRTY_LINES
  OUTPUT_STRIP_TRAILING_WHITESPACE)
