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
// Source:           ac_fir_const_coeffs_main_tb.cpp                                *
// Description:      C++ test bench                                                 *
//                                                                                  *
//***********************************************************************************

// To compile and run:
// $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include rtest_ac_fir_const_coeffs.cpp -o design
// ./design

#include <ac_dsp/ac_fir_const_coeffs.h>
#include <ac_channel.h>

#include <iostream>
using namespace std;
#include <fstream>

// Number of taps. Make sure that this number is the same as the number of
// elements in the wrapper class coeffs array.
const unsigned TAPS = 29;

const int  IN_W = 16;
const int  IN_I = 8;
const bool IN_S = true;

const int  OUT_W = 64;
const int  OUT_I = 32;
const bool OUT_S = true;

const int  COEFF_W = 32;
const int  COEFF_I = 16;
const bool COEFF_S = true;

const int  MAC_W = 64;
const int  MAC_I = 32;
const bool MAC_S = true;

const double PI_VALUE = 3.14159265358979323846;

// Declare types for input, output, coefficients and internal MAC variable
typedef ac_fixed < IN_W, IN_I, IN_S, AC_TRN, AC_WRAP > IN_TYPE_TB;
typedef ac_fixed < OUT_W, OUT_I, OUT_S, AC_TRN, AC_WRAP > OUT_TYPE_TB;
typedef ac_fixed < COEFF_W, COEFF_I, COEFF_S, AC_TRN, AC_WRAP > COEFF_TYPE_TB;
typedef ac_fixed < MAC_W, MAC_I, MAC_S, AC_TRN, AC_WRAP > MAC_TYPE_TB;

double abs_double(double input)
{
  return input < 0 ? -input : input;
}

//************************************************************************************************************************************************
//
// Class : ac_fir_const_coeffs_wrapper
// Description:
//   This class acts as a wrapper class around the ac_fir_const_coeffs class. Its main function is to hold the coeffs values as an array and
//   pass a pointer to the coeffs array in the constructor of the ac_fir_const_coeffs class, which is also a base class for the wrapper.
//   The user can configure the ac_fir_const_coeffs class to utilize different coefficient values by in turn utilizing differing wrapper classes
//   for every set of coefficient values they would like to test. Since this class is included in synthesis as the top level design, it ensures that
//   the ac_fir_const_coeffs testbench passes the synthesis and SCVerify stages.
//
//************************************************************************************************************************************************

#pragma hls_design top
template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, unsigned N_TAPS, FTYPE ftype >
class ac_fir_const_coeffs_wrapper : public ac_fir_const_coeffs<IN_TYPE, OUT_TYPE, COEFF_TYPE, ACC_TYPE, N_TAPS, ftype>
{
public:
  // The N_TAPS template parameter must be the same as the number of elements in the coeffs array below.
  // In this example, we are including the coeffs array from a text file, which is written in MATLAB. This is done
  // in order to make interfacing with matlab easier. The user can just as well copy-paste or manually type in the 
  // coefficient values themselves.
  const COEFF_TYPE coeffs[N_TAPS] = {
#include "ac_fir_const_coeffs_cfg.txt"
  };

  // Pass a pointer to the const coeffs array as a parameter in the ac_fir_const_coeffs constructor.
  ac_fir_const_coeffs_wrapper() : ac_fir_const_coeffs<IN_TYPE, OUT_TYPE, COEFF_TYPE, ACC_TYPE, N_TAPS, ftype> (coeffs) {  }
};

#define CW_CHECK_STREAM(istr, fname, desc) if (istr.fail()) { std::cerr << "Could not open " << desc << " file '" << fname << std::endl; return(-1); }

int main(int argc, char *argv[])
{
  fstream refout_fs;
  // Open text file that contains reference output data as an input file.
  refout_fs.open("ac_fir_const_coeffs_ref.txt", fstream::in);

  ac_channel < IN_TYPE_TB > input;
  ac_channel < OUT_TYPE_TB > output_test;
  OUT_TYPE_TB output_reference;

  // Enter frequencies of components to be tested as well as sampling frequency, in Hertz.
  // In order to avoid causing a problem when MATLAB tries to parse this cpp file, the user
  // is advised to leave the format of these declarations in their current form, only replacing the values with their own.
  // In other words, change the following line only if you know the ramifications of your changes.
  double F1 = 25, F2 = 150, Fs = 500;
  // Enter number of inputs.
  // This line, too, should only be changed if you know what you're doing, in order to ensure that MATLAB parses it well.
  // The user is advised to only change the value, not the way it is declared.
  const int n_inputs = 1024;
  // _input allows us to pass input data one-by-one through an ac_channel. _input_array contains
  // all the input data packed into an array.
  IN_TYPE_TB _input;
  IN_TYPE_TB _input_array[n_inputs];
  double input_type_max = _input.template set_val<AC_VAL_MAX>().to_double();
  double input_abs_max = 0.0;
  double mix_tone_inputs[n_inputs];

  for (int i = 0; i < n_inputs; i++) {
    // Mix sine waves with frequencies F1 and F2 to produce the mixed-tone input.
    mix_tone_inputs[i] = sin(2*PI_VALUE*F1*i/Fs) + sin(2*PI_VALUE*F2*i/Fs);
    // Store the maximum absolute value in the inputs, for use in normalization later.
    if (abs_double(mix_tone_inputs[i]) > input_abs_max) { input_abs_max = abs_double(mix_tone_inputs[i]); }
  }

  for (int i = 0; i < n_inputs; i++) {
    // Dividing every input signal produces a value that is normalized to 1.
    // Multiplying the division result by the maximum value representable for the input configuration
    // normalizes the input to that maximum value.
    _input_array[i] = (mix_tone_inputs[i] / input_abs_max) * input_type_max;
  }

  // Pass inputs to the top-level design.
  for (int i = 0; i < n_inputs; i++) {
    _input = _input_array[i];
    input.write(_input);
  }

  // Create a wrapper class object.
  ac_fir_const_coeffs_wrapper< IN_TYPE_TB, OUT_TYPE_TB, COEFF_TYPE_TB, MAC_TYPE_TB, TAPS, FOLD_ODD > filter;
  filter.run(input, output_test);

  double quant_noise_pow = 0, sig_pow = 0;
  double n_sample = 0;
  double output_test_design, output_test_ref;

  // Read data from output channels and find the sum of squares of signal and quantization noise values.
  while (output_test.available(1)) {
    n_sample++;
    output_test_design = output_test.read().to_double();
    refout_fs >> output_test_ref;
    // Use a typecast to convert the reference output to the output ac_fixed type, and then convert the
    // typecasted value back to double
    output_test_ref = ((OUT_TYPE_TB)output_test_ref).to_double();
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

  // If SQNR is above the pre-set threshold, the test passes. If it falls below the threshold, the test fails.
  if (SQNR < SQNR_TH) {
    printf("SQNR below threshold. Test FAILED.\n");
    return(-1);
  }

  printf("SQNR above threshold. Test PASSED.\n");
  return(0);
}

