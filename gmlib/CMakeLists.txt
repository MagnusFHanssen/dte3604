# ###############################################################################
# #
# # Copyright (C) 1994 Narvik University College
# # Contact: GMlib Online Portal at http://episteme.hin.no
# #
# # This file is part of the Geometric Modeling Library, GMlib.
# #
# # GMlib is free software: you can redistribute it and/or modify
# # it under the terms of the GNU Lesser General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.
# #
# # GMlib is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# # GNU Lesser General Public License for more details.
# #
# # You should have received a copy of the GNU Lesser General Public License
# # along with GMlib. If not, see <http://www.gnu.org/licenses/>.
# #
# ###############################################################################

# Minimum CMake version required
# -- 3.8  for compile featurs controll
# -- 3.11 for target declarations with no souces and a following
#         target_sources directive
cmake_minimum_required(VERSION 3.11)

# Nice name for the IDEs
project(gmlib VERSION 0.7.0 LANGUAGES CXX)

# IF on Wintendo platform
if(MSVC AND NOT VCPKG_TOOLCHAIN)
  message(FATAL_ERROR
    "Enable VCPKG toolchain on MSVC platform\n"
    "  Star over; CMAKE_TOOLCHAIN_FILE definition missing!\n"
    "  i.e.: cmake -DCMAKE_TOOLCHAIN_FILE=<path_to>\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake <path_to>\\source\\"
    )
endif()

# Add library target
if(MSVC)
  add_library(${PROJECT_NAME} STATIC)
else()
  add_library(${PROJECT_NAME} SHARED)
endif()


####################
# Configure Wintendo

# IF on Wintendo platform
if(MSVC AND NOT VCPKG_TOOLCHAIN)
  message(FATAL_ERROR
    "Enable VCPKG toolchain on MSVC platform\n"
    "  Star over; CMAKE_TOOLCHAIN_FILE definition missing!\n"
    "  i.e.: cmake -DCMAKE_TOOLCHAIN_FILE=<path_to>\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake <path_to>\\source\\"
    )
endif()



###########################
# Compiler spesific options

# Turn on c++11 compile features -- minimum CMake Version 3.8 required
target_compile_features(${PROJECT_NAME}
  INTERFACE $<$<CXX_COMPILER_ID:Clang>:cxx_std_11>
  INTERFACE $<$<CXX_COMPILER_ID:GNU>:cxx_std_11>
  INTERFACE $<$<CXX_COMPILER_ID:MSVC>:cxx_std_11>
  )

# Compile definitions
target_compile_definitions( ${PROJECT_NAME}
  INTERFACE $<$<CXX_COMPILER_ID:MSVC>: _USE_MATH_DEFINES>
  )

# Turn on/off internal opengl "safe"-call strategy
option(GMLIB_ENABLE_OPENGL_DEBUG_MODE "Enable opengl debug. OpenGL function calls are safe called." ON)
if(GMLIB_ENABLE_OPENGL_DEBUG_MODE)
  target_compile_definitions(${PROJECT_NAME}  PUBLIC GM_GL_DEBUG)
endif(GMLIB_ENABLE_OPENGL_DEBUG_MODE)

# Turn on/off stream operators in GMlib
option(GMLIB_ENABLE_STREAM_OPERATORS "Enable output stream operators for core types." ON)
if(GMLIB_ENABLE_STREAM_OPERATORS)
  target_compile_definitions(${PROJECT_NAME}  PUBLIC GM_STREAM)
endif(GMLIB_ENABLE_STREAM_OPERATORS)


# Compiler spesific options
target_compile_options(${PROJECT_NAME}
#  INTERFACE $<$<CXX_COMPILER_ID:AppleClang>:
#    -some-compiler-flag # somewhere over the rainbow
#    >
  INTERFACE $<$<CXX_COMPILER_ID:Clang>:
    -Weverything -pedantic -Werror
    -Wno-c++98-compat -Wno-c++98-compat-pedantic
    -Wno-documentation
    -Wno-exit-time-destructors
    -Wno-conversion
    -Wno-padded
    -Wno-weak-vtables
    -Wno-global-constructors
    -Wno-undefined-var-template
    -Wno-double-promotion
    -Wno-inconsistent-missing-override
    -Wno-deprecated
    -Wno-float-equal
    -Wno-ignored-qualifiers
    -Wno-non-virtual-dtor
    -Wno-covered-switch-default
    -Wno-gnu-statement-expression
    -Wno-conditional-uninitialized
    >
  INTERFACE $<$<CXX_COMPILER_ID:GNU>:
    -pedantic -Wall -Werror
    >
  INTERFACE $<$<CXX_COMPILER_ID:MSVC>:
    /w
    >
    )







target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/modules>
  )

set(SRCS_PREFIX modules)
set(CORE_SRCS_PREFIX ${SRCS_PREFIX}/core)
set(OPENGL_SRCS_PREFIX ${SRCS_PREFIX}/opengl)
set(SCENE_SRCS_PREFIX ${SRCS_PREFIX}/scene)
set(PARAMETRICS_SRCS_PREFIX ${SRCS_PREFIX}/parametrics)


# Core
target_sources(
  ${PROJECT_NAME} PRIVATE

  ${CORE_SRCS_PREFIX}/utils/gmcolor.cpp
  ${CORE_SRCS_PREFIX}/utils/gmstream.cpp
  )

# OpenGL
target_sources(
  ${PROJECT_NAME} PRIVATE

  ${OPENGL_SRCS_PREFIX}/gmbufferobject.cpp
  ${OPENGL_SRCS_PREFIX}/gmframebufferobject.cpp
  ${OPENGL_SRCS_PREFIX}/gmopenglmanager.cpp
  ${OPENGL_SRCS_PREFIX}/gmprogrampipeline.cpp
  ${OPENGL_SRCS_PREFIX}/gmprogram.cpp
  ${OPENGL_SRCS_PREFIX}/gmrenderbufferobject.cpp
  ${OPENGL_SRCS_PREFIX}/gmshader.cpp
  ${OPENGL_SRCS_PREFIX}/gmtexture.cpp

  ${OPENGL_SRCS_PREFIX}/bufferobjects/gmindexbufferobject.cpp
  ${OPENGL_SRCS_PREFIX}/bufferobjects/gmtexturebufferobject.cpp
  ${OPENGL_SRCS_PREFIX}/bufferobjects/gmuniformbufferobject.cpp
  ${OPENGL_SRCS_PREFIX}/bufferobjects/gmvertexbufferobject.cpp

  ${OPENGL_SRCS_PREFIX}/shaders/gmcomputeshader.cpp
  ${OPENGL_SRCS_PREFIX}/shaders/gmfragmentshader.cpp
  ${OPENGL_SRCS_PREFIX}/shaders/gmgeometryshader.cpp
  ${OPENGL_SRCS_PREFIX}/shaders/gmtesscontrolshader.cpp
  ${OPENGL_SRCS_PREFIX}/shaders/gmtessevaluationshader.cpp
  ${OPENGL_SRCS_PREFIX}/shaders/gmvertexshader.cpp
  )


# Scene
target_sources(
  ${PROJECT_NAME} PRIVATE

  ${SCENE_SRCS_PREFIX}/gmscene.cpp
  ${SCENE_SRCS_PREFIX}/gmsceneobject.cpp
  ${SCENE_SRCS_PREFIX}/gmvisualizer.cpp

  ${SCENE_SRCS_PREFIX}/camera/gmcamera.cpp
  ${SCENE_SRCS_PREFIX}/camera/gmisocamera.cpp

  ${SCENE_SRCS_PREFIX}/sceneobjects/gmcoordinatesystem.cpp
  ${SCENE_SRCS_PREFIX}/sceneobjects/gmpathtrack.cpp
  ${SCENE_SRCS_PREFIX}/sceneobjects/gmpathtrackarrows.cpp
  ${SCENE_SRCS_PREFIX}/sceneobjects/gmpointlightg.cpp
  ${SCENE_SRCS_PREFIX}/sceneobjects/gmspotlightg.cpp

  ${SCENE_SRCS_PREFIX}/event/gmevent.cpp
  ${SCENE_SRCS_PREFIX}/event/gmeventcontroller.cpp
  ${SCENE_SRCS_PREFIX}/event/gmeventmanager.cpp
  ${SCENE_SRCS_PREFIX}/event/gmsceneobjectattribute.cpp

  ${SCENE_SRCS_PREFIX}/light/gmlight.cpp
  ${SCENE_SRCS_PREFIX}/light/gmpointlight.cpp
  ${SCENE_SRCS_PREFIX}/light/gmspotlight.cpp
  ${SCENE_SRCS_PREFIX}/light/gmsun.cpp

  ${SCENE_SRCS_PREFIX}/render/gmdefaultrenderer.cpp
  ${SCENE_SRCS_PREFIX}/render/gmdefaultselectrenderer.cpp
  ${SCENE_SRCS_PREFIX}/render/gmrenderer.cpp
  ${SCENE_SRCS_PREFIX}/render/rendertargets/gmtexturerendertarget.cpp

  ${SCENE_SRCS_PREFIX}/utils/gmmaterial.cpp
  ${SCENE_SRCS_PREFIX}/utils/gmtexture.cpp

  ${SCENE_SRCS_PREFIX}/visualizers/gmcameravisualizer.cpp
  ${SCENE_SRCS_PREFIX}/visualizers/gmcoordsysrepvisualizer.cpp
  ${SCENE_SRCS_PREFIX}/visualizers/gmselectorvisualizer.cpp
  ${SCENE_SRCS_PREFIX}/visualizers/gmsurroundingspherevisualizer.cpp
  ${SCENE_SRCS_PREFIX}/visualizers/gmstdrepvisualizer.cpp
  ${SCENE_SRCS_PREFIX}/visualizers/gmpathtrackvisualizer.cpp
  ${SCENE_SRCS_PREFIX}/visualizers/gmpathtrackarrowsvisualizer.cpp
  )


# Parametrics -- no souces -- fully parameterized
#target_sources(
#  ${PROJECT_NAME} PRIVATE
#  )


# Trianglesystem -- no souces -- fully parameterized
#target_sources(
#  ${PROJECT_NAME} PRIVATE
#  )



##################
# Configure GLEW
find_package(GLEW REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC GLEW::GLEW )


##################
# Configure opengl
find_package(OpenGL REQUIRED)
if(MSVC)
  target_link_libraries(${PROJECT_NAME} PUBLIC GLU32 OpenGL32)
else()
  target_link_libraries(${PROJECT_NAME} PUBLIC GLU OpenGL)
endif()


##################
# Export targets and configuration
export(TARGETS ${PROJECT_NAME}
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}-targets.cmake"
  NAMESPACE gmlib::
  )

# Add alias gmlib::${PROJECT_NAME} to ${PROJECT_NAME} such that we can reference this
# in the same way as for other consumer usages (external library)
add_library(gmlib::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
    "cmake/${PROJECT_NAME}-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
    )



#########
# Testing
option(BUILD_TESTS "Build unit tests..." FALSE)
if(BUILD_TESTS)
  add_subdirectory(tests)
endif(BUILD_TESTS)



##############
# Benchmarking
option(BUILD_BENCHMARKS "Build benchmarks..." FALSE)
if(BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif(BUILD_BENCHMARKS)


