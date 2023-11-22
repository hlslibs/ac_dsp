/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.5                                                 *
 *                                                                        *
 *  Release Date    : Mon Nov 13 17:26:13 PST 2023                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.5.0                                               *
 *                                                                        *
 *  Copyright 2018 Siemens                                                *
 *                                                                        *
 **************************************************************************
 *  Licensed under the Apache License, Version 2.0 (the "License");       *
 *  you may not use this file except in compliance with the License.      * 
 *  You may obtain a copy of the License at                               *
 *                                                                        *
 *      http://www.apache.org/licenses/LICENSE-2.0                        *
 *                                                                        *
 *  Unless required by applicable law or agreed to in writing, software   * 
 *  distributed under the License is distributed on an "AS IS" BASIS,     * 
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or       *
 *  implied.                                                              * 
 *  See the License for the specific language governing permissions and   * 
 *  limitations under the License.                                        *
 **************************************************************************
 *                                                                        *
 *  The most recent version of this package is available at github.       *
 *                                                                        *
 *************************************************************************/
//***********************************************************************************
// Source:           ac_fir_prog_coeffs_main_tb.cpp                                 *
// Description:      C++ test bench                                                 *
//                                                                                  *
//***********************************************************************************

// To compile and run:
// $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_fir_prog_coeffs.cpp -o design
// ./design

#include <ac_dsp/ac_fir_prog_coeffs.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

const int TAPS = 27;
#define TYPE FOLD_ODD

// Declare types for input, output, coefficients and internal MAC variable
typedef ac_fixed < 28,  6, true, AC_TRN, AC_WRAP > IN_TYPE;
typedef ac_fixed < 64, 32, true, AC_TRN, AC_WRAP > OUT_TYPE;
typedef ac_fixed < 23,  7, true, AC_TRN, AC_WRAP > COEFF_TYPE;
typedef ac_fixed < 64, 32, true, AC_TRN, AC_WRAP > MAC_TYPE;

const double PI_VALUE = 3.14159265358979323846;

using namespace std;
#define CW_CHECK_STREAM(istr, fname, desc) if (istr.fail()) { std::cerr << "Could not open " << desc << " file '" << fname << std::endl; return(-1); }

int main(int argc, char *argv[])
{
  fstream refout_fs, coeff_fs;

  // Open text files that contain reference output and coefficents data as input files.
  refout_fs.open("ac_fir_prog_coeffs_ref.txt", fstream::in);
  coeff_fs.open("ac_fir_prog_coeffs_cfg.txt", fstream::in);

  ac_channel < IN_TYPE > input;
  ac_channel < OUT_TYPE > output_test;
  ac_fir_prog_coeffs < IN_TYPE, OUT_TYPE, COEFF_TYPE, MAC_TYPE, TAPS, TYPE > filter;
  ac_channel < COEFF_TYPE > coeffs_ch;
  COEFF_TYPE coeffs[TAPS];
  OUT_TYPE output_reference;
  // Enter frequencies of components to be tested as well as sampling frequency, in Hertz.
  double F1 = 25, F2 = 150, Fs = 500;
  // Enter number of inputs.
  const int n_inputs = 1024;
  // _input allows us to pass input data one-by-one through an ac_channel. _input_array contains
  // all the input data packed into an array.
  IN_TYPE _input, _input_array[n_inputs];
  double input_type_max = _input.template set_val<AC_VAL_MAX>().to_double();
  double input_abs_max = 0.0;
  double mix_tone_inputs[n_inputs];

  for (int i = 0; i < n_inputs; i++) {
    // Mix sine waves with frequencies F1 and F2 to produce the mixed-tone input.
    mix_tone_inputs[i] = sin(2*PI_VALUE*F1*((double)i)/Fs) + sin(2*PI_VALUE*F2*((double)i)/Fs);
    // Store the maximum absolute value in the inputs, for use in normalization later.
    if (abs(mix_tone_inputs[i]) > input_abs_max) { input_abs_max = abs(mix_tone_inputs[i]); }
  }

  for (int i = 0; i < n_inputs; i++) {
    // Dividing every input signal produces a value that is normalized to 1.
    // Multiplying the division result by the maximum value representable for the input configuration
    // normalizes the input to that maximum value.
    _input_array[i] = (mix_tone_inputs[i] / input_abs_max) * input_type_max;
  }

  double _coeff;

  // Read coefficients data from the corresponding text file and load to top-level design.
  for (int i = 0; i < TAPS; i++) {
    coeff_fs >> _coeff;
    coeffs[i] = _coeff;
  }

  // Pass inputs to the top-level design.
  for (int i = 0; i < n_inputs; i++) {
    _input = _input_array[i];
    input.write(_input);
    filter.run(input, output_test, coeffs);
  }

  double quant_noise_pow = 0, sig_pow = 0;
  double n_sample = 0;
  double output_test_design, output_test_ref;

  // Read data from output channels and find the sum of squares of signal and quantization noise values.
  while (output_test.available(1)) {
    n_sample++;
    output_test_design = output_test.read().to_double();
    refout_fs >> output_test_ref;

    quant_noise_pow += (output_test_design - output_test_ref)*(output_test_design - output_test_ref);
    sig_pow += output_test_ref*output_test_ref;
  }

  // Define an SQNR threshold over here, in dB.
  double SQNR_TH = 60;
  // Note that there is no need to take RMS voltage into consideration while calculating SQNR, the
  // denominator (i.e. number of signal samples) gets canceled out in the formula below. Just dividing
  // the sum of squares of signal values by the sum of squares of quantization noise values should be enough.
  double SQNR = 10*log10(sig_pow/quant_noise_pow);

  cout << "SQNR = " << SQNR << "dB" << endl;

  refout_fs.close();
  coeff_fs.close();

  // If SQNR is above the pre-set threshold, the test passes. If it falls below the threshold, the test fails.
  if (SQNR < SQNR_TH) {
    printf("SQNR below threshold. Test FAILED.\n");
    return(-1);
  }

  printf("SQNR above threshold. Test PASSED.\n");
  return(0);
}
