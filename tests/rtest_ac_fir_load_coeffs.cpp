//***********************************************************************************
// Source:           ac_fir_load_coeffs_main_tb.cppp                                *
// Description:      C++ test bench                                                 *
//                                                                                  *
//***********************************************************************************

// To compile and run:
// $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include -I$MGC_HOME/shared/pkgs/ac_dsp/include/ac_dsp ac_fir_load_coeffs_main_tb.cpp -o design
// ./design

#include <ac_fixed.h>
#include <ac_channel.h>
#include <ac_dsp/ac_fir_load_coeffs.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

// Number of taps
const unsigned TAPS = 27;

const int  IN_W = 32;
const int  IN_I = 16;
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

// Declare types for input, output, coefficients and internal MAC variable
typedef ac_fixed < IN_W, IN_I, IN_S, AC_TRN, AC_WRAP > IN_TYPE_TB;
typedef ac_fixed < OUT_W, OUT_I, OUT_S, AC_TRN, AC_WRAP > OUT_TYPE_TB;
typedef ac_fixed < COEFF_W, COEFF_I, COEFF_S, AC_TRN, AC_WRAP > COEFF_TYPE_TB;
typedef ac_fixed < MAC_W, MAC_I, MAC_S, AC_TRN, AC_WRAP > MAC_TYPE_TB;

using namespace std;
#define CW_CHECK_STREAM(istr, fname, desc) if (istr.fail()) { std::cerr << "Could not open " << desc << " file '" << fname << std::endl; return(-1); }

int main(int argc, char *argv[])
{
  fstream refout_fs, coeff_fs;

  // Open text files that contain reference output and coefficents data as input files.
  refout_fs.open("ac_fir_load_coeffs_ref.txt", fstream::in);
  coeff_fs.open("ac_fir_load_coeffs_cfg.txt", fstream::in);

  ac_channel<IN_TYPE_TB> input;
  ac_channel<OUT_TYPE_TB> output_test;
  ac_channel<COEFF_TYPE_TB> coeffs_ch;
  ac_channel<bool> load;

  COEFF_TYPE_TB coeffs_arr_elem;
  OUT_TYPE_TB output_reference;
  // Enter frequencies of components to be tested as well as sampling frequency, in Hertz.
  // In order to avoid causing a problem when the MATLAB script tries to parse this cpp file, the user
  // is advised to leave the format of these declarations in their current form, only replacing the values with their own.
  // In other words, change the following line only if you know the ramifications of your changes.
  double F1 = 25, F2 = 150, Fs = 500;
  // Enter number of inputs.
  // This line, too, should only be changed if you know what you're doing, in order to ensure that MATLAB parses it well.
  // The user is advised to only change the value, not the way it is declared.
  const int n_inputs = 1024;
  // _input allows us to pass input data one-by-one through an ac_channel. _input_array contains
  // all the input data packed into an array.
  IN_TYPE_TB _input, _input_array[n_inputs];
  double input_type_max = _input.template set_val<AC_VAL_MAX>().to_double();
  double input_abs_max = 0.0;
  double mix_tone_inputs[n_inputs];

  for (int i = 0; i < n_inputs; i++) {
    // Mix sine waves with frequencies F1 and F2 to produce the mixed-tone input.
    mix_tone_inputs[i] = sin(2*M_PI*F1*i/Fs) + sin(2*M_PI*F2*i/Fs);
    // Store the maximum absolute value in the inputs, for use in normalization later.
    if (abs(mix_tone_inputs[i]) > input_abs_max) { input_abs_max = abs(mix_tone_inputs[i]); }
  }

  for (int i = 0; i < n_inputs; i++) {
    // Dividing every input signal by the max value in the signal produces a value that is normalized to 1.
    // Multiplying the division result by the maximum value representable for the input bitwidths
    // normalizes the input to that maximum value.
    _input_array[i] = (mix_tone_inputs[i] / input_abs_max) * input_type_max;
  }

  double _coeff;

  // Load coefficients data into the top level design, coefficient-by-coefficient
  // Read coefficients data from the corresponding text file and load to top-level design.
  for (int i = 0; i < TAPS; i++) {
    coeff_fs >> _coeff;
    coeffs_arr_elem = _coeff;
    coeffs_ch.write(coeffs_arr_elem);
  }

  // First load coefficients into the design.
  ac_fir_load_coeffs < IN_TYPE_TB, OUT_TYPE_TB, COEFF_TYPE_TB, MAC_TYPE_TB, TAPS, FOLD_ODD > filter;
  load.write(true);
  filter.run(input, coeffs_ch, output_test, load);
  load.write(false);

  // Now that the coefficients are loaded, pass inputs to the top-level design.
  for (int i = 0; i < n_inputs; i++) {
    _input = _input_array[i];
    input.write(_input);
  }

  // Run top level filter function again in order to process all inputs.
  filter.run(input, coeffs_ch, output_test, load);

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
  coeff_fs.close();

  // If SQNR is above the pre-set threshold, the test passes. If it falls below the threshold, the test fails.
  if (SQNR < SQNR_TH) {
    printf("SQNR below threshold. Test FAILED.\n");
    return(-1);
  }

  printf("SQNR above threshold. Test PASSED.\n");
  return(0);
}
