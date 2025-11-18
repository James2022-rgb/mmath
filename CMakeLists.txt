cmake_minimum_required(VERSION 3.22)

# --------------------------------------------------------------------------------
# Internal variables
#

set(TARGET_NAME "mmath")
set(SRC_DIR "src")

# --------------------------------------------------------------------------------
# Project
#

project(${TARGET_NAME}
  VERSION 0.1.0.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=thread-safety")
endif()

add_library(${TARGET_NAME} STATIC)
target_compile_features(${TARGET_NAME} PRIVATE cxx_std_23)

if(MSVC)
  target_compile_options(${TARGET_NAME} PRIVATE /W4)
else()
  target_compile_options(${TARGET_NAME} PRIVATE -Wall -Wextra -Werror)
endif()

#
#
#

set(HEADERS_MMATH_PUBLIC_DIR "${SRC_DIR}/mmath/public")
set(HEADERS_MMATH_PUBLIC
  ${HEADERS_MMATH_PUBLIC_DIR}/geometry.h
  ${HEADERS_MMATH_PUBLIC_DIR}/matrix.h
  ${HEADERS_MMATH_PUBLIC_DIR}/projection.h
  ${HEADERS_MMATH_PUBLIC_DIR}/quat.h
  ${HEADERS_MMATH_PUBLIC_DIR}/vector.h
)
source_group("Public" FILES ${HEADERS_MMATH_PUBLIC})

set(SOURCES_MMATH_PRIVATE_DIR "${SRC_DIR}/mmath/private")
set(SOURCES_MMATH_PRIVATE
  ${SOURCES_MMATH_PRIVATE_DIR}/geometry.cpp
  ${SOURCES_MMATH_PRIVATE_DIR}/matrix.cpp
  ${SOURCES_MMATH_PRIVATE_DIR}/projection.cpp
  ${SOURCES_MMATH_PRIVATE_DIR}/quat.cpp
  ${SOURCES_MMATH_PRIVATE_DIR}/vector.cpp
)
source_group("Private" FILES ${SOURCES_MMATH_PRIVATE})

set(SOURCES
  ${SOURCES_MMATH_PRIVATE}
)

set(HEADERS
  ${HEADERS_MMATH_PUBLIC}
)

target_sources(${TARGET_NAME} PRIVATE ${SOURCES})
target_sources(${TARGET_NAME} PRIVATE ${HEADERS})

target_include_directories(${TARGET_NAME} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${SOURCES_MMATH_PRIVATE_DIR}
)

target_include_directories(${TARGET_NAME} PUBLIC
  ${SRC_DIR}
)

target_link_libraries(${TARGET_NAME} PUBLIC
  mbase
)
