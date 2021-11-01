# Find GTest
find_package(GTest MODULE REQUIRED)

# Enable testing
enable_testing()

# Configure threads
find_package(Threads REQUIRED)

# Test sources
set(TESTS
  core_static_staticproc_compiletests
  core_containers_gmarray_tests
  core_containers_dvectorn_tests
  scene_sceneobject_tests
  parametrics_curves_compiletests
  parametrics_surfaces_compiletests
  parametrics_transform_tests
  parametrics_object_creation_tests
  )


# Add tests
foreach(TEST ${TESTS})

  add_executable( unittest_${TEST} ${TEST}.cc )

  target_link_libraries(unittest_${TEST}
    PRIVATE
      gmlib::gmlib
      GTest::GTest GTest::Main
      Threads::Threads
    )

  # Turn off platform-spesific extensions
  set_target_properties(unittest_${TEST} PROPERTIES CXX_EXTENSIONS OFF)

  # Add additional compile options
  target_compile_options(unittest_${TEST}
    PUBLIC $<$<CXX_COMPILER_ID:Clang>:
      -Wno-used-but-marked-unused
      >
      )

  gtest_add_tests( unittest_${TEST} "" AUTO )
endforeach(TEST)

message("GMlib unittests SETUP done!")
