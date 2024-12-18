/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) Math Library                                       *
 *                                                                        *
 *  Software Version: 3.6                                                 *
 *                                                                        *
 *  Release Date    : Tue Nov 12 23:14:00 PST 2024                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.6.0                                               *
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
// Source:           rtest_ac_cic_dec_full.cpp                                      *
// Description:      C++ test bench                                                 *
//                                                                                  *
//***********************************************************************************

// To compile and run:
//   $MGC_HOME/bin/g++ -I$MGC_HOME/shared/include rtest_ac_cic_dec_full.cpp -o ac_cic_dec_full
//   ./ac_cic_dec_full

#include <ac_dsp/ac_cic_dec_full.h>
#include <ac_channel.h>
#include <fstream>
#include <ac_fixed.h>

#include "ac_cic_dec_full_param.h"

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
  unsigned err_cnt = 0;
  fstream inf, ref;
  cout << "=============================================================================" << endl;
  cout << "------------------------------ Running  rtest_ac_cic_dec_full.cpp ----------------------------------" << endl;
  cout << "=============================================================================" << endl;
  cout << "Testing class: ac_cic_dec_full" << endl;

  string input_file = "ac_cic_dec_full_input.txt";
  string output_ref = "ac_cic_dec_full_ref.txt";

  ac_channel < IN_TYPE_TB > fixed_in;
  ac_channel < OUT_TYPE_TB > fixed_out;
  ac_channel < long double > fixed_ref;

  double big_diff = 0, data_in, diff;
  long double ref_in;

  // Matlab pads a zero at the beginning of the output reference vector. In order
  // to compensate for this behavior of Matlab, we pass in an extra zero at the beginning
  // of the input vector to the design. The counters and flags that are initialized within the design
  // library files themselves are initialized in such a way that the first data element passed
  // at the input is passed as-is, straight to the output. Hence the extra zero we pass to
  // the input over here is directly reflected at the design output as well.
  data_in = 0;
  fixed_in.write(data_in);

  // Read input values
  cout << "Reading input vectors from: " << input_file << endl;
  inf.open(input_file.c_str(), fstream::in);
  CW_CHECK_STREAM(inf,input_file,"Input data");
  do {
    inf >> data_in;
    if (! inf.eof()) fixed_in.write(data_in);
  } while (! inf.eof());
  inf.close();

  // Read reference values
  cout << "Reading output reference from: " << output_ref << endl;
  ref.open(output_ref.c_str(), fstream::in);
  CW_CHECK_STREAM(ref,output_ref,"Matlab reference output");
  do {
    long double ref_in;
    ref >> ref_in;
    if (!ref.eof()) fixed_ref.write(ref_in);
  } while (! ref.eof());
  ref.close();

  cout << "Number of inputs = " << fixed_in.debug_size() << endl;
  cout << "Number of reference outputs = " << fixed_ref.debug_size() << endl;

  ac_cic_dec_full < IN_TYPE_TB, OUT_TYPE_TB, R_TB, M_TB, N_TB > filter;
  filter.run(fixed_in, fixed_out);

  cout << "Number of DUT outputs = " << fixed_ref.debug_size() << endl << endl;

  // Make sure that the design outputs more values than/the same number of values as the reference.
  if (fixed_out.debug_size() < fixed_ref.debug_size()) { 
    cout << "Fewer outputs produced : reference = " << fixed_ref.debug_size() << " values, DUT = " << fixed_out.debug_size() << " values." << endl;
    err_cnt++;
  }

  double n_output = 0;

  while (fixed_out.available(1) && fixed_ref.available(1)) {
    n_output++;
    OUT_TYPE_TB fixed_data_out = fixed_out.read(); // get DUT output

    long double ref_in = fixed_ref.read(); // get reference output
    diff = custom_abs(ref_in - fixed_data_out.to_double());

    if (diff != 0) {
      cout << "	@ERROR: Data mismatch on sample " << n_output << " : expected   " << ref_in << "  Received : " << fixed_data_out << endl;
      err_cnt++;
    }

    if (diff > big_diff) {
      big_diff = diff;
    }
  }

  if (err_cnt) {
    cout << "Test FAILED. Worst case absolute difference = " << big_diff << endl;
  } else {
    cout << "Test PASSED." << endl;
  }

  return(err_cnt);
}
