/************************************************************************************/
/*                                 CatWare                                          */
/*                 Copyright (c) 2003-2016 Mentor Graphics Corporation.             */
/*                             All Rights Reserved                                  */
/*                                                                                  */
/* This document contains information that is proprietary to Mentor Graphics        */
/* Corporation. The original recipient of this document may duplicate this          */
/* document in whole or in part for internal business purposes only, provided       */
/* that this entire notice appears in all copies. In duplicating any part of        */
/* this document, the recipient agrees to make every reasonable effort to           */
/* prevent the unauthorized use and distribution of the proprietary information.    */
/*                                                                                  */
/************************************************************************************/
/************************************************************************************/
/* NO WARRANTY. MENTOR GRAPHICS CORPORATION EXPRESSLY DISCLAIMS ALL WARRANTY        */
/* FOR THE SOFTWARE. TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE                  */
/* LAW, THE SOFTWARE AND ANY RELATED DOCUMENTATION IS PROVIDED "AS IS"              */
/* AND WITH ALL FAULTS AND WITHOUT WARRANTIES OR CONDITIONS OF ANY                  */
/* KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, THE              */
/* IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR                  */
/* PURPOSE, OR NONINFRINGEMENT. THE ENTIRE RISK ARISING OUT OF USE OR               */
/* DISTRIBUTION OF THE SOFTWARE REMAINS WITH YOU.                                   */
/*                                                                                  */
/************************************************************************************/
/* THE SOFTWARE IS PROVIDED WITH NO SUPPORT OR INDEMNIFICATION.                     */
/************************************************************************************/

/************************************************************************************/
/* Source:           ac_cic_dec_full_main_tb.cpp                                    *
 * Description:      C++ test bench                                                 *
 *                                                                                  */
/************************************************************************************/

// To compile and run:
// $MGC_HOME/bin/c++ -std=c++11 -I$MGC_HOME/shared/include -I$MGC_HOME/shared/pkgs/ac_dsp/include/ac_dsp ac_cic_dec_full_main_tb.cpp -o design
// ./design

#include <ac_dsp/ac_cic_dec_full.h>
#include <ac_channel.h>
#include <fstream>
#include <ac_fixed.h>

const unsigned R_TB = 7;
const unsigned M_TB = 2;
const unsigned N_TB = 4;

const int  in_W = 32;
const int  in_I = 16;
const bool in_S = true;

const int  out_W = 48;
const int  out_I = 32;
const bool out_S = true;

typedef ac_fixed <in_W, in_I, in_S> IN_TYPE_TB;
typedef ac_fixed <out_W, out_I, out_S> OUT_TYPE_TB;

#define NO_IN_DATA 10000

#include <iostream>
using namespace std;

#define CW_CHECK_STREAM(istr, fname, desc) if (istr.fail()) { cerr << "Could not open " << desc << " file '" << fname << endl; return(-1); }

// Define a custom absolute value function for double values. This is done because using the built-in "abs" function for doubles sometimes causes naming conflicts.
double custom_abs(double input)
{
  return input < 0 ? -input : input;
}

int main(int argc, char *argv[])
{
  fstream inf, ref;
  cout << "=============================================================================" << endl;
  cout << "------------------------------ Running  rtest_ac_cic_dec_full.cpp ----------------------------------" << endl;
  cout << "=============================================================================" << endl;
  cout << "Testing class: ac_cic_dec_full" << endl;

  string input_file = "ac_cic_dec_full_input.txt";
  string output_ref = "ac_cic_dec_full_ref.txt";

  cout << "Reading input vectors from: " << input_file << endl;
  cout << "Reading output reference from: " << output_ref << endl;

  inf.open(input_file.c_str(), fstream::in);   /* Channel 0 Input data */
  ref.open(output_ref.c_str(), fstream::in);    /* Channel 0 reference output */
  CW_CHECK_STREAM(inf,input_file,"Channel 0 input data");
  CW_CHECK_STREAM(ref,output_ref,"Channel 0 reference output");

  ac_channel < IN_TYPE_TB > fixed_in;
  ac_channel < OUT_TYPE_TB > fixed_out;

  double big_diff = 0, data_in, diff;
  long double ref_in;

  data_in = 0;
  fixed_in.write(data_in);

  for (int inread = 0; inread < NO_IN_DATA; inread++) {
    inf >> data_in;
    fixed_in.write(data_in);
  }

  ac_cic_dec_full < IN_TYPE_TB, OUT_TYPE_TB, R_TB, M_TB, N_TB > filter;
  filter.run(fixed_in, fixed_out);

  bool test_pass = true;
  double n_output = 0;

  while (fixed_out.available(1)) {
    n_output++;
    OUT_TYPE_TB fixed_data_out = fixed_out.read();

    ref >> ref_in;
    diff = custom_abs(ref_in - fixed_data_out.to_double());

    if (diff != 0) {
#ifdef DEBUG
      cout << "	@ERROR: Data mismatch : expected   " << ref_in << "  Received : " << fixed_data_out << endl;
#endif
      test_pass = false;
    }

    if (diff > big_diff) {
      big_diff = diff;
    }
  }

  cout << "Number of inputs  = " << NO_IN_DATA << endl;
  cout << "Number of outputs = " << n_output << endl << endl;

  inf.close();
  ref.close();

  if (!test_pass) {
    cout << "Test FAILED. Worst case absolute difference = " << big_diff << endl;
    return(-1);
  }

  cout << "Test PASSED." << endl;
  return(0);

}
