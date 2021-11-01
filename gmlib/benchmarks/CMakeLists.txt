# Find benchmark
find_package(benchmark REQUIRED)

# Enable testing
enable_testing()

# Configure threads
find_package(Threads REQUIRED)

# Test sources
set(BENCHMARKS
  core_containers_array_benchmarks
  parametrics_pcurve_evaluate_benchmarks
  )

# Add tests
foreach(BENCH ${BENCHMARKS})

  add_executable( bench_${BENCH} ${BENCH}.cc )

  target_link_libraries(bench_${BENCH}
    PRIVATE ${PROJECT_NAME}
      gmlib::gmlib
      benchmark::benchmark benchmark::benchmark_main
      Threads::Threads
      )

  # Turn off platform-spesific extensions
  set_target_properties(bench_${BENCH} PROPERTIES CXX_EXTENSIONS OFF)

  # Add additional compile options
#    target_compile_options(bench_${BENCH}
#      PUBLIC $<$<CXX_COMPILER_ID:Clang>:
#        -Wno-some-compile-feature
#        >
#        )

  add_test(NAME bench_${BENCH} COMMAND bench_${BENCH})
endforeach(BENCH)

message("GMlib benchmarks SETUP done!")


