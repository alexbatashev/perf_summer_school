#include "task2.hpp"

#include <benchmark/benchmark.h>
#include <random>
#include <vector>

inline constexpr int TEST_SIZE = 50000;

static void evt_dispatch(benchmark::State &state) {
  std::random_device rd{};
  auto eng = std::mt19937{rd()};
  auto gen = std::uniform_int_distribution{0, 63};

  std::vector<EventID> ids;
  ids.reserve(100);

  EventDispatcher d;

  for (auto _ : state) {
    if (!d.hasFreeSlots()) {
      size_t idx = gen(eng);
      d.release(ids[idx]);
      ids.erase(std::next(ids.begin(), idx));
    }
    
    ids.push_back(d.acquire());
  }
}

// Register the function as a benchmark
BENCHMARK(evt_dispatch)->Unit(benchmark::kMicrosecond);

// Run the benchmark
BENCHMARK_MAIN();
