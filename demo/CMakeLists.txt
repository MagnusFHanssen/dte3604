# Minimum CMake version required
# -- 3.8  for compile featurs controll
# -- 3.11 for target declarations with no souces and a following
#         target_sources directive
cmake_minimum_required(VERSION 3.11)

# Global name of the project
project(QmlDemo VERSION 0.2.0 LANGUAGES CXX)

# Add executable target
add_executable( ${PROJECT_NAME} )



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


####################
# Configure gmlib2qt
find_package(gmlib 0.7 REQUIRED CONFIG NO_DEFAULT_PATH)
target_link_libraries( ${PROJECT_NAME} gmlib::gmlib )


################################
# Configure Qt5 packages
find_package(Qt5Core 5.1)
find_package(Qt5Quick 5.1)
find_package(Qt5Gui 5.1)
target_link_libraries(  ${PROJECT_NAME}
  Qt5::Core Qt5::Quick Qt5::Gui
  )


###########################
# Compiler spesific options

# Get GMlib Transitive target COMPILE properties
get_target_property(DEMO_GMLIB_TRANS_INTERFACE_COMPILE_DEFINITIONS gmlib::gmlib INTERFACE_COMPILE_DEFINITIONS)
get_target_property(DEMO_GMLIB_TRANS_INTERFACE_COMPILE_FEATURES gmlib::gmlib INTERFACE_COMPILE_FEATURES)
get_target_property(DEMO_GMLIB_TRANS_INTERFACE_COMPILE_OPTIONS gmlib::gmlib INTERFACE_COMPILE_OPTIONS)

# Set GMlib INTERFACE target COMPILE property dependency on GMlib2 PUBLIC profile
target_compile_definitions(${PROJECT_NAME} PUBLIC ${DEMO_GMLIB_TRANS_INTERFACE_COMPILE_DEFINITIONS})
target_compile_features(${PROJECT_NAME}    PUBLIC ${DEMO_GMLIB_TRANS_INTERFACE_COMPILE_FEATURES})
target_compile_options(${PROJECT_NAME}     PUBLIC ${DEMO_GMLIB_TRANS_INTERFACE_COMPILE_OPTIONS})

# Turn off platform-spesific extensions
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

# Add additional compile options
target_compile_options(${PROJECT_NAME}
#  PUBLIC $<$<CXX_COMPILER_ID:AppleClang>:
#    -some-compiler-flag # somewhere over the rainbow
#    >
  PUBLIC $<$<CXX_COMPILER_ID:Clang>:
    -Wno-unused-parameter
    -Wno-old-style-cast
    -Wno-implicit-fallthrough
    -Wno-shadow
    -Wno-reserved-id-macro
    -Wno-unreachable-code
    -Wno-unreachable-code-return
    -Wno-unused-macros
    -Wno-undefined-reinterpret-cast
    -Wno-redundant-parens
    -Wno-inconsistent-missing-destructor-override
    -Wno-zero-as-null-pointer-constant
    -Wno-extra-semi-stmt
    >
#  PUBLIC $<$<CXX_COMPILER_ID:GNU>:
#    -some-compiler-flag
#    >
#  PUBLIC $<$<CXX_COMPILER_ID:MSVC>:
#    /some-compiler-option
#    >
    )




# Qt moc'ing
qt5_wrap_cpp( HDRS_MOC
  hidmanager/hidaction.h
  hidmanager/hidmanager.h
  hidmanager/hidmanagertreemodel.h
  hidmanager/standardhidmanager.h
  hidmanager/defaulthidmanager.h

  application/fboinsgrenderer.h
  application/gmlibwrapper.h
  application/guiapplication.h
  application/window.h

  scenario.h
  )
target_sources( ${PROJECT_NAME} PRIVATE ${HDRS_MOC})



# Qt resouces files
qt5_add_resources( RCCS_MOC
  application/application.qrc
  )
target_sources( ${PROJECT_NAME} PRIVATE ${RCCS_MOC})


# Sources
target_sources( ${PROJECT_NAME} PRIVATE
  hidmanager/hidaction.cpp
  hidmanager/hidbinding.cpp
  hidmanager/hidinput.cpp
  hidmanager/hidinputevent.cpp
  hidmanager/hidkbmouseinput.cpp
  hidmanager/hidmanager.cpp
  hidmanager/hidmanagertreemodel.cpp
  hidmanager/standardhidmanager.cpp
  hidmanager/defaulthidmanager.cpp

  application/fboinsgrenderer.cpp
  application/gmlibwrapper.cpp
  application/guiapplication.cpp
  application/window.cpp

  application/main.cpp

  scenario.cpp
  )

