// This example is borrowed from perf-ninja course
// https://github.com/dendibakh/perf-ninja/blob/main/labs/core_bound/compiler_intrinsics_1/solution.cpp

#include "task1.hpp"

#ifdef __riscv
#define VECTORIZED __attribute__((target("arch=rv64gcv")))
#elif defined(__x86_64)
#define VECTORIZED __attribute__((target("arch=skylake")))
#else
#define VECTORIZED
#endif

#include <cmath>
#include <cstdint>
#include <random>

constexpr int N = 40000;

using int8x8_t = uint8_t __attribute__((ext_vector_type(8)));
using int16x8_t = uint16_t __attribute__((ext_vector_type(8)));

template<typename T>
__attribute__((always_inline)) constexpr T shift_left(T x, int len) {
  T result = 0;
  const int vecLen = __builtin_vectorelements(T);

  for (int i = 0; i < vecLen - len; i++) {
    result[i+len] = x[i];
  }

  return result;
}

template<typename T>
constexpr T scan_add(T x) {
  int vecLen = __builtin_vectorelements(T);
  int numIter = __builtin_ctz(vecLen);

  T shift = x;
  int mul = 1;

  #pragma loop unroll(full)
  for (int i = 0; i < numIter; i++) {
      shift = shift + shift_left(shift, mul);
      mul <<= 1;
  }

  return shift;
}

template<typename T, typename P>
T vload(const P *ptr) {
  int vecLen = __builtin_vectorelements(T);

  T result{0};
  for (int i = 0; i < vecLen; i++) {
    result[i] = ptr[i];
  }

  return result;
}

template<typename T, typename P> void vstore(P *ptr, T x) {
  int vecLen = __builtin_vectorelements(T);

  for (int i = 0; i < vecLen; i++) {
    ptr[i] = x[i];
  }
}

VECTORIZED void imageSmoothing_vector(const InputVector &input, uint8_t radius,
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

  constexpr int VL = 8;
  const uint8_t* subtractPtr = input.data() + pos - radius - 1;
  const uint8_t* addPtr = input.data() + pos + radius;
  uint16_t* outputPtr = output.data() + pos;
  int16x8_t current = currentSum;

  int i = 0;
  for (; i + 7 < limit - pos; i += VL) {
    int8x8_t sub_8 = vload<int8x8_t>(subtractPtr + i);
    int16x8_t sub = __builtin_convertvector(sub_8, int16x8_t);
    int8x8_t add_8 = vload<int8x8_t>(addPtr + i);
    int16x8_t add = __builtin_convertvector(add_8, int16x8_t);

    int16x8_t diff = add - sub;

    int16x8_t scan = scan_add(diff);
    int16x8_t result = scan + current;

    vstore(outputPtr + i, result);

    currentSum = result[7];
    current = currentSum;
  }

  pos += i;

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

void imageSmoothing_generic(const InputVector &input, uint8_t radius,
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

static void *imageSmoothing_resolver() {
  return (void*)&imageSmoothing_vector;
  // return (void*)&imageSmoothing_generic;
}

__attribute__((ifunc("_ZL23imageSmoothing_resolverv")))
void imageSmoothing(const InputVector &input, uint8_t radius,
                    OutputVector &output);

void zero(OutputVector &data, std::size_t size) {
  data.clear();
  data.resize(size);
}
