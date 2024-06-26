#include "task2.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>

inline constexpr int TEST_SIZE = 100;

using Catch::Matchers::WithinRel;

inline constexpr std::array<double, TEST_SIZE * 2 - 1> GOLDEN = {
    -1,       -1.96,    -2.8804,  -3.7616,  -4.604,   -5.408,   -6.174,
    -6.9024,  -7.5936,  -8.248,   -8.866,   -9.448,   -9.9944,  -10.5056,
    -10.982,  -11.424,  -11.832,  -12.2064, -12.5476, -12.856,  -13.132,
    -13.376,  -13.5884, -13.7696, -13.92,   -14.04,   -14.13,   -14.1904,
    -14.2216, -14.224,  -14.198,  -14.144,  -14.0624, -13.9536, -13.818,
    -13.656,  -13.468,  -13.2544, -13.0156, -12.752,  -12.464,  -12.152,
    -11.8164, -11.4576, -11.076,  -10.672,  -10.246,  -9.7984,  -9.3296,
    -8.84,    -8.33,    -7.8,     -7.2504,  -6.6816,  -6.094,   -5.488,
    -4.864,   -4.2224,  -3.5636,  -2.888,   -2.196,   -1.488,   -0.7644,
    -0.0256,  0.728,    1.496,    2.278,    3.0736,   3.8824,   4.704,
    5.538,    6.384,    7.2416,   8.1104,   8.99,     9.88,     10.78,
    11.6896,  12.6084,  13.536,   14.472,   15.416,   16.3676,  17.3264,
    18.292,   19.264,   20.242,   21.2256,  22.2144,  23.208,   24.206,
    25.208,   26.2136,  27.2224,  28.234,   29.248,   30.264,   31.2816,
    32.3004,  33.32,    32.34,    31.36,    30.3804,  29.4016,  28.424,
    27.448,   26.474,   25.5024,  24.5336,  23.568,   22.606,   21.648,
    20.6944,  19.7456,  18.802,   17.864,   16.932,   16.0064,  15.0876,
    14.176,   13.272,   12.376,   11.4884,  10.6096,  9.74,     8.88,
    8.03,     7.1904,   6.3616,   5.544,    4.738,    3.944,    3.1624,
    2.3936,   1.638,    0.896,    0.168,    -0.5456,  -1.2444,  -1.928,
    -2.596,   -3.248,   -3.8836,  -4.5024,  -5.104,   -5.688,   -6.254,
    -6.8016,  -7.3304,  -7.84,    -8.33,    -8.8,     -9.2496,  -9.6784,
    -10.086,  -10.472,  -10.836,  -11.1776, -11.4964, -11.792,  -12.064,
    -12.312,  -12.5356, -12.7344, -12.908,  -13.056,  -13.178,  -13.2736,
    -13.3424, -13.384,  -13.398,  -13.384,  -13.3416, -13.2704, -13.17,
    -13.04,   -12.88,   -12.6896, -12.4684, -12.216,  -11.932,  -11.616,
    -11.2676, -10.8864, -10.472,  -10.024,  -9.542,   -9.0256,  -8.4744,
    -7.888,   -7.266,   -6.608,   -5.9136,  -5.1824,  -4.414,   -3.608,
    -2.764,   -1.8816,  -0.9604};

TEST_CASE("fft", "[task2]") {
  // generate random polynomials
  polynomial<double> poly1(TEST_SIZE);
  polynomial<double> poly2(TEST_SIZE);
  polynomial<double> poly3(TEST_SIZE * 2 - 1);

  for (int n = 0; n < TEST_SIZE; ++n) {
    poly1[n] = -1.0 + 0.02 * n;
    poly2[n] = 1.0 - 0.02 * n;
  }

  poly3 = poly1 * poly2;

  for (int n = 0; n < TEST_SIZE * 2 - 1; ++n) {
    REQUIRE_THAT(poly3[n], WithinRel(GOLDEN[n], 0.01));
  }
}