#include <benchmark/benchmark.h>

#include <parametrics/curves/gmpcircle.h>
using namespace GMlib;



static void BM_PCircle_Resample(benchmark::State& state)
{
  auto pcircle = PCircle<float>();
  DVector<DVector<Vector<float,3>>> samps;

  // The test loop
  while (state.KeepRunning()) {
    pcircle.resample(samps, state.range(0), state.range(1));
  }
}


BENCHMARK(BM_PCircle_Resample)
  ->Unit(benchmark::kNanosecond)
  ->RangeMultiplier(2)
  ->Ranges({{2, 2 << 15}, {1,2}});
