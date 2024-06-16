#include "task2.hpp"

#include <benchmark/benchmark.h>
#include <random>

inline constexpr int TEST_SIZE = 50000;

static void fft(benchmark::State &state) {
  std::random_device rd{};
  auto eng = std::mt19937{rd()};
  auto gen = std::uniform_real_distribution{-1.0, 1.0};

  // generate random polynomials
  polynomial<double> poly1(TEST_SIZE);
  polynomial<double> poly2(TEST_SIZE);
  polynomial<double> poly3(TEST_SIZE * 2 - 1);

  for (int n = 0; n < TEST_SIZE; ++n) {
    poly1[n] = gen(eng);
    poly2[n] = gen(eng);
  }

  for (auto _ : state) {
    poly3 = poly1 * poly2;
    benchmark::DoNotOptimize(poly3);
  }
}

// Register the function as a benchmark
BENCHMARK(fft)->Unit(benchmark::kMicrosecond);

// Run the benchmark
BENCHMARK_MAIN();