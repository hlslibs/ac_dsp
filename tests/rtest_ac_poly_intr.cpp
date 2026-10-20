/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <ac_dsp/ac_poly_intr.h>
#include <iostream>

using namespace std;

typedef ac_fixed < 8, 4, true, AC_TRN, AC_WRAP > IN_TYPE;
typedef ac_fixed < 8, 2, true > COEFF_TYPE;
typedef ac_fixed < 24, 12, true > OUT_TYPE;

struct STR_CTRL_TYPE {
  bool sign[3];
  ac_int < 8, false > corr[3];
};

template < int COEFFSZ >
struct STR_COEFF_TYPE {
  COEFF_TYPE coeffs[COEFFSZ];
};

template < int NTAPS, int COEFFSZ, FTYPE ftype >
int test_impulses(const COEFF_TYPE (&coeffs)[COEFFSZ], const int (&fir)[3 * NTAPS], const char *name)
{
  int errors = 0;
  for (int raw = -128; raw < 128; raw++) {
    ac_poly_intr < IN_TYPE, COEFF_TYPE, OUT_TYPE, OUT_TYPE, STR_CTRL_TYPE,
                  STR_COEFF_TYPE < COEFFSZ >, NTAPS, COEFFSZ, 3, ftype > filter;
    ac_channel < IN_TYPE > input;
    ac_channel < OUT_TYPE > output;
    ac_channel < STR_CTRL_TYPE > ctrl;
    ac_channel < STR_COEFF_TYPE < COEFFSZ > > coeff;
    ac_channel < bool > load;

    STR_CTRL_TYPE c = {{true, true, false}, {2, 1, 0}};
    STR_COEFF_TYPE < COEFFSZ > h;
    for (int i = 0; i < COEFFSZ; i++) {
      h.coeffs[i] = coeffs[i];
    }
    ctrl.write(c);
    coeff.write(h);
    load.write(true);
    filter.run(input, output, ctrl, coeff, load);

    IN_TYPE impulse;
    impulse.set_slc(0, ac_int < 8, true > (raw));
    // One extra zero input drains the one-input latency.
    for (int n = 0; n <= NTAPS; n++) {
      input.write(n == 0 ? impulse : IN_TYPE(0));
      load.write(false);
      filter.run(input, output, ctrl, coeff, load);
    }

    bool passed = true;
    for (int n = 0; n < 3 * NTAPS; n++) {
      if (!output.available(1)) {
        cout << name << ": missing output " << n << " for input " << raw << endl;
        passed = false;
        break;
      }
      // Original FIR coefficients are in eighths, input codes in sixteenths,
      // and output codes in units of 2^-12. This is an exact impulse oracle.
      int expected = raw * fir[n] * 32;
      int actual = (output.read() * 4096).to_int();
      if (actual != expected) {
        cout << name << ": input " << raw << ", output " << n
             << ": expected " << expected << ", received " << actual << endl;
        passed = false;
      }
    }
    if (output.available(1)) {
      cout << name << ": unexpected extra output for input " << raw << endl;
      passed = false;
    }
    if (!passed) { errors++; }
  }
  cout << name << ": 256 input codes, " << errors << " failed." << endl;
  return errors;
}

int main()
{
  const COEFF_TYPE even_coeffs[3] = {0.5, 0.25, -0.5};
  const COEFF_TYPE odd_coeffs[6] = {0.5, 0.25, 0.25, 0.5, -0.5, 0};
  // Original FIR responses, before the symmetric-pair coefficient transform.
  const int even_fir[6] = {4, 2, 0, 0, 2, 4};
  const int odd_fir[9] = {4, 2, 0, 1, 4, 1, 0, 2, 4};

  int errors = test_impulses < 2, 3, FOLD_EVEN > (even_coeffs, even_fir, "FOLD_EVEN");
  errors += test_impulses < 3, 6, FOLD_ODD > (odd_coeffs, odd_fir, "FOLD_ODD");
  cout << (errors ? "Test FAILED." : "Test PASSED.") << endl;
  return errors ? -1 : 0;
}
