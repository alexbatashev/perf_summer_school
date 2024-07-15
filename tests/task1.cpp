#include "task1.hpp"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

static void reference_solution(const InputVector &input, uint8_t radius,
                               OutputVector &output) {
  int pos = 0;
  int currentSum = 0;
  int size = (int)input.size();

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

  // 2. main loop. During optimization, focus mainly on this part
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

TEST_CASE("Smoke test", "[task1]") {
  InputVector inA;
  init(inA);

  OutputVector expected, received;
  zero(expected, (int)inA.size());
  zero(received, (int)inA.size());

  reference_solution(inA, radius, expected);
  imageSmoothing(inA, radius, received);

  REQUIRE(expected.size() == received.size());

  // std::cout << "expected: ";
  // for (auto v : expected)
  //   std::cout << v << " ";
  // std::cout << "\n";
  //
  // std::cout << "received: ";
  // for (auto v : received)
  //   std::cout << v << " ";
  // std::cout << "\n";

  auto cmp_result =
      std::mismatch(expected.begin(), expected.end(), received.begin());
  REQUIRE(cmp_result.first == expected.end());
}
