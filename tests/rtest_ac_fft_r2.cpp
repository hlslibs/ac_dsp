// Licensed under the Apache License, Version 2.0.
// See LICENSE for details.

// Compile with -DTEST_DIF and/or -DTEST_SDF to select the FFT architecture.
#if defined(TEST_DIF) && defined(TEST_SDF)
#include <ac_dsp/ac_fft_dif_r2_sdf.h>
#elif defined(TEST_DIF)
#include <ac_dsp/ac_fft_dif_r2_inpl.h>
#elif defined(TEST_SDF)
#include <ac_dsp/ac_fft_dit_r2_sdf.h>
#else
#include <ac_dsp/ac_fft_dit_r2_inpl.h>
#endif

#include <iostream>
using namespace std;

template <int TWID_PREC>
bool test_driver(const int expected[8][2])
{
  typedef ac_fixed < 32, 4, true > dType;
  typedef ac_complex < dType > cx_dType;
#if defined(TEST_DIF) && defined(TEST_SDF)
  ac_fft_dif_r2_sdf < 8, 32, TWID_PREC, 32, 4 > fft;
#elif defined(TEST_DIF)
  ac_fft_dif_r2_inpl < 8, TWID_PREC, 32, 4 > fft;
#elif defined(TEST_SDF)
  ac_fft_dit_r2_sdf < 8, 32, TWID_PREC, 32, 4 > fft;
#else
  ac_fft_dit_r2_inpl < 8, TWID_PREC, 32, 4 > fft;
#endif
  ac_channel < cx_dType > input, output;

  // Natural input x[1] = 1. DIT requires bit-reversed input.
  for (int i = 0; i < 8; i++) {
    cx_dType x(0, 0);
#ifdef TEST_DIF
    if (i == 1) { x.r() = 1; }
#else
    if (i == 4) { x.r() = 1; }
#endif
    input.write(x);
  }
  fft.run(input, output);
#ifdef TEST_SDF
  // Drain the one-frame SDF startup latency.
  for (int i = 0; i < 8; i++) { input.write(cx_dType(0, 0)); }
  fft.run(input, output);
#endif

  if (output.size() != 8) {
    cout << "TWID_PREC = " << TWID_PREC << ": expected 8 outputs, got " << output.size() << endl;
    return false;
  }

  bool passed = true;
  for (int i = 0; i < 8; i++) {
#ifdef TEST_DIF
    // DIF returns bit-reversed output.
    int bin = ((i & 1) << 2) | (i & 2) | ((i & 4) >> 2);
#else
    int bin = i;
#endif
    cx_dType y = output.read();
    int real = y.r().template slc<32>(0).to_int();
    int imag = y.i().template slc<32>(0).to_int();
    if (real != expected[bin][0] || imag != expected[bin][1]) {
      cout << "TWID_PREC = " << TWID_PREC << ", bin " << bin
           << ": expected (" << expected[bin][0] << ", " << expected[bin][1]
           << "), got (" << real << ", " << imag << ")" << endl;
      passed = false;
    }
  }
  cout << "TWID_PREC = " << TWID_PREC << ": Test " << (passed ? "PASSED." : "FAILED.") << endl;
  return passed;
}

int main()
{
  // Exact Q4.28 outputs from integer arithmetic with each stage's rounding
  // and cumulative 1/8 scaling. Nearest-rounded sqrt(2)/2 twiddle raw values
  // are 92682, 185364, 2965821, 189812531 and 759250125 for these widths.
  const int expected_20[8][2] = {
    {33554432, 0}, {23726592, -23726592},
    {0, -33554432}, {-23726592, -23726592},
    {-33554432, 0}, {-23726592, 23726592},
    {0, 33554432}, {23726592, 23726592},
  };
  const int expected_24[8][2] = {
    {33554432, 0}, {23726568, -23726568},
    {0, -33554432}, {-23726568, -23726568},
    {-33554432, 0}, {-23726568, 23726568},
    {0, 33554432}, {23726568, 23726568},
  };
#ifdef TEST_DIF
  const int expected_30[8][2] = {
    {33554432, 0}, {23726566, -23726567},
    {0, -33554432}, {-23726566, -23726567},
    {-33554432, 0}, {-23726566, 23726567},
    {0, 33554432}, {23726566, 23726567},
  };
  const int expected_32[8][2] = {
    {33554432, 0}, {23726566, -23726567},
    {0, -33554432}, {-23726567, -23726567},
    {-33554432, 0}, {-23726566, 23726567},
    {0, 33554432}, {23726567, 23726567},
  };
#else
  const int expected_30[8][2] = {
    {33554432, 0}, {23726566, -23726567},
    {0, -33554432}, {-23726567, -23726567},
    {-33554432, 0}, {-23726566, 23726566},
    {0, 33554432}, {23726566, 23726566},
  };
#endif
  bool passed = true;
  passed = test_driver<19>(expected_20) && passed;
  passed = test_driver<20>(expected_20) && passed;
  passed = test_driver<24>(expected_24) && passed;
  passed = test_driver<30>(expected_30) && passed;
#ifdef TEST_DIF
  passed = test_driver<32>(expected_32) && passed;
#else
  passed = test_driver<32>(expected_30) && passed;
#endif
  return passed ? 0 : 1;
}
