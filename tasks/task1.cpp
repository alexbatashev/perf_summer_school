// This example is borrowed from perf-ninja course
// https://github.com/dendibakh/perf-ninja/blob/main/labs/core_bound/compiler_intrinsics_1/solution.cpp

#include "task1.hpp"

#include <cmath>
#include <limits>
#include <random>

constexpr int N = 40000;

void imageSmoothing(const InputVector &input, uint8_t radius,
                    OutputVector &output) {
  int pos = 0;
  int currentSum = 0;
  int size = static_cast<int>(input.size());

  // 1. left border - time spend in this loop can be ignored, no need to
  // optimize it
  for (int i = 0; i < std::min<int>(size, radius); ++i) {
    currentSum += input[i];
  }

  int limit = std::min(radius + 1, size - radius);
  for (pos = 0; pos < limit; ++pos) {
    currentSum += input[pos + radius];
    output[pos] = currentSum;
  }

  // 2. main loop.
  limit = size - radius;
  for (; pos < limit; ++pos) {
    currentSum -= input[pos - radius - 1];
    currentSum += input[pos + radius];
    output[pos] = currentSum;
  }

  // 3. special case, executed only if size <= 2*radius + 1
  limit = std::min(radius + 1, size);
  for (; pos < limit; pos++) {
    output[pos] = currentSum;
  }

  // 4. right border - time spend in this loop can be ignored, no need to
  // optimize it
  for (; pos < size; ++pos) {
    currentSum -= input[pos - radius - 1];
    output[pos] = currentSum;
  }
}

void init(InputVector &data) {
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, 255);

  data.reserve(N);
  for (int i = 0; i < N; i++) {
    uint8_t value = static_cast<uint8_t>(distribution(generator));
    data.emplace_back(value);
  }
}

void zero(OutputVector &data, std::size_t size) {
  data.clear();
  data.resize(size);
}