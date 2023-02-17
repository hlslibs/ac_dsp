/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) DSP Library                                        *
 *                                                                        *
 *  Software Version: 3.4                                                 *
 *                                                                        *
 *  Release Date    : Mon Feb  6 09:12:03 PST 2023                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.4.6                                               *
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
//***********************************************************************************************************//
// File:         ac_poly_intr.h
// Description:
//    This filter is designed using class based hierarchy. The member function "run()" contained in the
//    ac_poly_intr class is the top function for the C++ design.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_poly_intr.h>
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Define input, output, coefficient and accumulator type
//      typedef ac_fixed<32, 16, true> IN_TYPE;
//      typedef ac_fixed<64, 32, true> OUT_TYPE;
//      typedef ac_fixed<32, 16, true> COEFF_TYPE;
//      typedef ac_fixed<64, 32, true> ACC_TYPE;
//      typedef struct {
//        COEFF_TYPE coeffs[NTAPS * DF];
//      }
//      typedef struct {
//      bool sign[IR];
//      ac_int < 8, false > corr[IR];
//      }
//      STR_CTRL_TYPE;
//      STR_COEFF_TYPE;
//      ac_poly_intr < IN_TYPE, COEFF_TYPE, ACC_TYPE, OUT_TYPE, STR_CTRL_TYPE, STR_COEFF_TYPE, NTAPS, COEFFSZ,
//      IR, TYPE > filter_design1;
//      ac_channel < IN_TYPE > fixed_in;
//      ac_channel < OUT_TYPE > fixed_out;
//      ac_channel < bool > sign_st;
//      ac_channel < STR_CTRL_TYPE > ctrl_st;
//      ac_channel < STR_COEFF_TYPE > coeffs_st;
//      STR_CTRL_TYPE  ctrl_st_t;
//      STR_COEFF_TYPE coeffs_st_t;
//      // Make sure to initialize the inputs and the coeffs
//      // Call the top-level function.
//      filter_design1.run(input, output, ctrl_st, coeffs_st);
//
//      CCS_RETURN(0);
//    }
//
// Revision History:
//    3.4.0  - Added CDesignChecker waivers/fixes for ac_dsp IP blocks.
//             Changes made in general:
//               - CNS violations were waived away.
//               - FXD violations were fixed by assigning fixed value
//                 initializations to ac_fixed variables.
//               - ABR and UMR violations were waived away.
//
//***********************************************************************************************************//

#ifndef _INCLUDED_AC_POLY_INTR_H_
#define _INCLUDED_AC_POLY_INTR_H_

#include <ac_fixed.h>
#include <ac_int.h>
#include <ac_channel.h>

// Make sure that the enum is only defined once.
#ifndef __POLY_FILTER_TYPES_ENUM_DEF__
#define __POLY_FILTER_TYPES_ENUM_DEF__
// The parameters within this enum help the user choose between different filter architectures.
typedef enum { FOLD_EVEN, FOLD_ODD, FOLD_ANTI } FTYPE;

#endif

//=================================================================================================================
// Class: ac_poly_intr_core
// Description: The class member functions implement different architectures for the polyphase interpolator filter.
//-----------------------------------------------------------------------------------------------------------------
template < class IN_TYPE, class COEFF_TYPE, class ACC_TYPE, class OUT_TYPE, int NTAPS, int COEFFSZ, int IF >
class ac_poly_intr_core
{
private: // Data
  IN_TYPE taps[NTAPS];
  ACC_TYPE acc_a[IF];
  ACC_TYPE acc_b[IF];
  bool flip, init;

public: // Functions
  // Constructor
  ac_poly_intr_core() {
    ac::init_array < AC_VAL_0 > (acc_a, IF);
    ac::init_array < AC_VAL_0 > (acc_b, IF);
    ac::init_array < AC_VAL_0 > (taps, NTAPS);
    flip = false;
    init = false;
  }

//------------------------------------------------------------------------------------------------------
// Member Function: ac_polyIntrSymmetricEvenTaps()
// Description: ac_polyIntrSymmetricEvenTaps() implements a symmetric filter with even number of Taps
//
  void ac_polyIntrSymmetricEvenTaps ( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, COEFF_TYPE coeffs[COEFFSZ], bool sign[IF], ac_int < 8, false > corr[IF]) {
#pragma hls_unroll
    SHIFT_REG: for (int i = (NTAPS - 1); i >= 1; i--) {
      taps[i] = taps[i - 1];
    }
    INTR_F: for (int j = 0; j < IF; j++) {
      if (j == 0) {
        taps[0] = data_in.read();
        flip = !flip;
      }
      ACC_TYPE fold, acc;
      acc = 0.0;
#pragma hls_unroll
      MAC_E: for (int i = (NTAPS / 2) - 1; i >= 0; i--) {
        IN_TYPE tp;
        if (sign[j]) {                                      // sign specifies whether filter is symmetric or anti-symmetric
          tp = taps[NTAPS - 1 - i];
        } else {
          tp = -taps[NTAPS - 1 - i];
        }
        fold = (taps[i] + tp);                             // fold the input for symmetry
        acc += coeffs[i + j * NTAPS / 2] * fold;
      }
      ACC_TYPE t1, t2;
      if (flip) {
        acc_b[j] = acc;
        t1 = acc_a[j];
        t2 = acc_a[corr[j]];
      } else {
        acc_a[j] = acc;
        t1 = acc_b[j];
        t2 = acc_b[corr[j]];
      }
      if (init)
        if (j != corr[j]) {                               // if coefficients are modified based on symmetric pair technique
          ACC_TYPE tn;
          if (sign[j]) {
            tn = -t2;
          } else {
            tn = t2;
          }
          data_out.write((t1 + tn) >> 1);
        } else {                                          // Coefficients are symmetric without modifications
          data_out.write(t1);
        }
    }
    init = true;
  }

//------------------------------------------------------------------------------------------------------
// Member Function: ac_polyIntrSymmetricOddTaps()
// Description: ac_polyIntrSymmetricOddTaps() implements a symmetric filter with even number of Taps
//
  void ac_polyIntrSymmetricOddTaps ( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, COEFF_TYPE coeffs[COEFFSZ], bool sign[IF], ac_int < 8, false > corr[IF]) {
#pragma hls_unroll
    SHIFT_REG: for (int i = (NTAPS - 1); i >= 1; i--) {
      taps[i] = taps[i - 1];
    }
    INTR_F: for (int j = 0; j < IF; j++) {
      if (j == 0) {
        taps[0] = data_in.read();
        flip = !flip;
      }
      ACC_TYPE fold, acc;
      acc = 0.0;
#pragma hls_unroll
      MAC_O: for (int i = 0; i < (((NTAPS - 1) / 2) + 1); i++) {
        if (i == (NTAPS - 1) / 2) {
          fold = taps[i]; // Central tap passes through for odd symmetric*/
        } else {
          IN_TYPE tp;
          if (sign[j]) {                                   // sign specifies whether filter is symmetric or anti-symmetric
            tp = taps[NTAPS - 1 - i];
          } else {
            tp = -taps[NTAPS - 1 - i];
          }
          fold = (taps[i] + tp);  // fold the input for symmetry
        }
        acc += coeffs[i + (NTAPS / 2 + 1) * j] * fold;
      }
      ACC_TYPE t1, t2;
      if (flip) {
        acc_b[j] = acc;
        t1 = acc_a[j];
#pragma hls_waive ABR
        t2 = acc_a[corr[j]];
      } else {
        acc_a[j] = acc;
        t1 = acc_b[j];
#pragma hls_waive ABR
        t2 = acc_b[corr[j]];
      }
      if (init)
        if (j != corr[j]) {                               // if coefficients are modified based on symmetric pair technique
          ACC_TYPE tn;
          if (sign[j]) {
            tn = -t2;
          } else {
            tn = t2;
          }
          data_out.write((t1 + tn) >> 1);
        } else {                                          // Coefficients are symmetric without modifications
          data_out.write(t1);
        }
    }                                                            // end of IF loop
    init = true;
  }

//------------------------------------------------------------------------------------------------------
// Member Function: ac_polyIntrAntiSymmetric()
// Description: ac_polyIntrAntiSymmetric() implements a symmetric filter with even number of Taps  
//
  void ac_polyIntrAntiSymmetric ( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, COEFF_TYPE coeffs[COEFFSZ], bool sign[IF], ac_int < 8, false > corr[IF]) {
#pragma hls_unroll
    SHIFT_REG: for (int i = (NTAPS - 1); i >= 1; i--) {
      taps[i] = taps[i - 1];
    }
    INTR_F: for (int j = 0; j < IF; j++) {
      if (j == 0) {
        taps[0] = data_in.read();
        flip = !flip;
      }
      ACC_TYPE fold, acc;
      acc = 0.0;
#pragma hls_unroll
      MAC_N: for (int i = (NTAPS - 1); i >= 0; i--) {
        acc += taps[i] * coeffs[i + NTAPS * j];
      }
      data_out.write(acc);
    }                                                            // end of IF loop
    init = true;
  }

};

// Include header for CCS_BLOCK
#include <mc_scverify.h>

//===================================================================================================================
// Class: ac_poly_intr
// Description: This class contains the top function for the polyphase interpolator filter ('run()'). 
//-------------------------------------------------------------------------------------------------------------------
template < class IN_TYPE, class COEFF_TYPE, class ACC_TYPE, class OUT_TYPE, class STR_CTRL_TYPE, class STR_COEFF_TYPE, int NTAPS, int COEFFSZ, int IF, FTYPE ftype >
class ac_poly_intr
{

private: // Data
  ac_poly_intr_core < IN_TYPE, COEFF_TYPE, ACC_TYPE, OUT_TYPE, NTAPS, COEFFSZ, IF > filter_inst;
  STR_CTRL_TYPE ctrl_t;
  STR_COEFF_TYPE coeffs_t;

public: // Functions
  // Constructor
  ac_poly_intr() {
    //read_control_chan = false;
  }

//------------------------------------------------------------------------------------------------------
// Member Function: run()
// Description: run() is top function for C++ module.  
//
#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void CCS_BLOCK(run)( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, ac_channel < STR_CTRL_TYPE > &ctrl_st, ac_channel < STR_COEFF_TYPE > &coeffs_st, ac_channel < bool > &read_ctrl_chan) {
		bool read_ctrl = read_ctrl_chan.read();
    if (read_ctrl) {
      ctrl_t = ctrl_st.read();
      coeffs_t = coeffs_st.read();
    } else {
#pragma hls_waive CNS
      if (ftype == FOLD_EVEN) {
#pragma hls_waive UMR
        filter_inst.ac_polyIntrSymmetricEvenTaps(data_in, data_out, coeffs_t.coeffs, ctrl_t.sign, ctrl_t.corr); // Symmetric filter with even number of Taps
      }
#pragma hls_waive CNS
      if (ftype == FOLD_ODD) {
#pragma hls_waive UMR
        filter_inst.ac_polyIntrSymmetricOddTaps(data_in, data_out, coeffs_t.coeffs, ctrl_t.sign, ctrl_t.corr);  // Symmetric filter with odd number of Taps
      }
#pragma hls_waive CNS
      if (ftype == FOLD_ANTI) {
#pragma hls_waive UMR
        filter_inst.ac_polyIntrAntiSymmetric(data_in, data_out, coeffs_t.coeffs, ctrl_t.sign, ctrl_t.corr);     // Antisymmetric filter
      }
    }
  }
};

#endif

