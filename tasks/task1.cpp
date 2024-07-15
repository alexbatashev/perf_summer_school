// This example is borrowed from perf-ninja course
// https://github.com/dendibakh/perf-ninja/blob/main/labs/core_bound/compiler_intrinsics_1/solution.cpp

#include "task1.hpp"

#include <iostream>
#include <utility>

#ifdef __riscv
#include <riscv_vector.h>
#define VECTORIZED
// #define VECTORIZED __attribute__((target("arch=rv64gcv")))
#elif defined(__x86_64)
#define VECTORIZED __attribute__((target("arch=skylake")))
#else
#define VECTORIZED
#endif

#include <cmath>
#include <cstdint>
#include <random>

constexpr int N = 40000;

namespace detail {
template <typename T, size_t Size>
struct simd_native {
#if defined(__clang__)
  using type = T __attribute__((ext_vector_type(Size)));
#endif
};

#if !defined(__clang__)
template <>
struct simd_native<uint8_t, 8> {
#ifdef __riscv
  using type = vuint8m1_t;
#else
  using type = uint8_t __attribute__((vector_size(8 * sizeof(uint8_t))));
#endif
};

template <>
struct simd_native<uint16_t, 8> {
#ifdef __riscv
  using type = vuint16m1_t;
#else
  using type = uint16_t __attribute__((vector_size(8 * sizeof(uint16_t))));
#endif
};
#endif
}

template <typename T>
concept scalar = std::integral<T> || std::floating_point<T>;

template <scalar T, size_t Size>
class simd {
public:
  using native_type = typename detail::simd_native<T, Size>::type;

  simd() = default;

  explicit simd(T scalar) : data_(native_type{scalar}) {}

  constexpr size_t size() const {
    return Size;
  }

  T operator[](size_t idx) const {
    return data_[idx];
  }

  template <typename U>
  simd<U, Size> convert() const {
    auto cvt = __builtin_convertvector(data_, typename detail::simd_native<U, Size>::type);
    return simd<U, Size>(cvt);
  }

  static simd load(const T *ptr) {
    simd<T, Size> result;

    #pragma GCC ivdep
    for (size_t i = 0; i < Size; i++) {
      result.data_[i] = ptr[i];
    }

    return result;
  }

  void store(T *ptr) const {
    #pragma GCC ivdep
    for (size_t i = 0; i < Size; i++) {
      ptr[i] = data_[i];
    }
  }

  simd<T, Size> shift_left(size_t len) const {
    native_type result{0};

    for (int i = 0; i < size() - len; i++) {
      result[i+len] = data_[i];
    }

    return simd<T, Size>(result);
  }

  simd<T, Size> operator-(const simd<T, Size> &other) {
    native_type res = data_ - other.data_;
    return simd<T, Size>(res);
  }

  simd<T, Size> operator+(const simd<T, Size> &other) {
    native_type res = data_ + other.data_;
    return simd<T, Size>(res);
  }

private:
  template <scalar U, size_t S>
  friend class simd;
  template <scalar U, size_t S, typename Op>
  simd<U, S> scan(const simd<T, S> x, Op op);

  simd(native_type data) : data_(data) {}

  native_type data_{0};
};

template <scalar T, size_t Size, typename Op = std::plus<>>
simd<T, Size> scan(const simd<T, Size> x, Op op = Op()) {
  size_t numIter = __builtin_ctz(x.size());

  auto shift = x;
  int mul = 1;

  #pragma clang loop unroll(full)
  for (int i = 0; i < numIter; i++) {
    shift = shift + shift.shift_left(mul);
    mul <<= 1;
  }

  return shift;
}

VECTORIZED void imageSmoothing_vector(const InputVector &input, uint8_t radius,
                    OutputVector &output) {
  int pos = 0;
  uint16_t currentSum = 0;
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

  constexpr size_t VL = 8;
  const uint8_t* subtractPtr = input.data() + pos - radius - 1;
  const uint8_t* addPtr = input.data() + pos + radius;
  uint16_t* outputPtr = output.data() + pos;
  simd<uint16_t, 8> current{currentSum};

  int i = 0;
  for (; i + 7 < limit - pos; i += VL) {
    auto sub_8 = simd<uint8_t, VL>::load(subtractPtr + i);
    auto sub = sub_8.convert<uint16_t>();
    auto add_8 = simd<uint8_t, VL>::load(addPtr + i);
    auto add = add_8.convert<uint16_t>();

    auto diff = add - sub;
    auto scan = ::scan(diff);

    auto result = scan + current;
    result.store(outputPtr + i);

    currentSum = result[7];
    current = simd<uint16_t, VL>(currentSum);
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
