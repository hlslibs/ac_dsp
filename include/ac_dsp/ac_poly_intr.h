////////////////////////////////////////////////////////////////////////////////
// Catapult Synthesis
// 
// Copyright (c) 2003-2018 Mentor Graphics Corp.
//       All Rights Reserved
// 
// This document contains information that is proprietary to Mentor Graphics
// Corp. The original recipient of this document may duplicate this  
// document in whole or in part for internal business purposes only, provided  
// that this entire notice appears in all copies. In duplicating any part of  
// this document, the recipient agrees to make every reasonable effort to  
// prevent the unauthorized use and distribution of the proprietary information.
//
////////////////////////////////////////////////////////////////////////////////

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

  void ac_polyIntrSymmetricEvenTaps ( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, COEFF_TYPE coeffs[COEFFSZ], bool sign[IF], ac_int < 8, false > corr[IF]) {
#pragma unroll
    SHIFT_REG: for (int i = (NTAPS - 1); i >= 1; i--) {
      taps[i] = taps[i - 1];
    }
    INTR_F: for (int j = 0; j < IF; j++) {
      if (j == 0) {
        taps[0] = data_in.read();
        flip = !flip;
      }
      ACC_TYPE fold, acc;
      acc = 0;
#pragma unroll
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

  void ac_polyIntrSymmetricOddTaps ( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, COEFF_TYPE coeffs[COEFFSZ], bool sign[IF], ac_int < 8, false > corr[IF]) {
#pragma unroll
    SHIFT_REG: for (int i = (NTAPS - 1); i >= 1; i--) {
      taps[i] = taps[i - 1];
    }
    INTR_F: for (int j = 0; j < IF; j++) {
      if (j == 0) {
        taps[0] = data_in.read();
        flip = !flip;
      }
      ACC_TYPE fold, acc;
      acc = 0;
#pragma unroll
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
    }                                                            // end of IF loop
    init = true;
  }

  void ac_polyIntrAntiSymmetric ( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, COEFF_TYPE coeffs[COEFFSZ], bool sign[IF], ac_int < 8, false > corr[IF]) {
#pragma unroll
    SHIFT_REG: for (int i = (NTAPS - 1); i >= 1; i--) {
      taps[i] = taps[i - 1];
    }
    INTR_F: for (int j = 0; j < IF; j++) {
      if (j == 0) {
        taps[0] = data_in.read();
        flip = !flip;
      }
      ACC_TYPE fold, acc;
      acc = 0;
#pragma unroll
      MAC_N: for (int i = (NTAPS - 1); i >= 0; i--) {
        acc += taps[i] * coeffs[i + NTAPS * j];
      }
      data_out.write(acc);
    }                                                            // end of IF loop
    init = true;
  }

};

template < class IN_TYPE, class COEFF_TYPE, class ACC_TYPE, class OUT_TYPE, class STR_CTRL_TYPE, class STR_COEFF_TYPE, int NTAPS, int COEFFSZ, int IF, FTYPE ftype >
class ac_poly_intr
{

private: // Data
  ac_poly_intr_core < IN_TYPE, COEFF_TYPE, ACC_TYPE, OUT_TYPE, NTAPS, COEFFSZ, IF > filter_inst;
  STR_CTRL_TYPE ctrl_t;
  STR_COEFF_TYPE coeffs_t;

public: // Functions
  // Constructor
  ac_poly_intr() {}

#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void run( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, ac_channel < STR_CTRL_TYPE > &ctrl_st, ac_channel < STR_COEFF_TYPE > &coeffs_st) {
#ifndef __SYNTHESIS__
    while (ctrl_st.available(1))
#endif
    {
      ctrl_t = ctrl_st.read();
    }

#ifndef __SYNTHESIS__
    while (coeffs_st.available(1))
#endif
      coeffs_t = coeffs_st.read();

    if (ftype == FOLD_EVEN) {
      filter_inst.ac_polyIntrSymmetricEvenTaps(data_in, data_out, coeffs_t.coeffs, ctrl_t.sign, ctrl_t.corr); // Symmetric filter with even number of Taps
    }
    if (ftype == FOLD_ODD) {
      filter_inst.ac_polyIntrSymmetricOddTaps(data_in, data_out, coeffs_t.coeffs, ctrl_t.sign, ctrl_t.corr);  // Symmetric filter with odd number of Taps
    }
    if (ftype == FOLD_ANTI) {
      filter_inst.ac_polyIntrAntiSymmetric(data_in, data_out, coeffs_t.coeffs, ctrl_t.sign, ctrl_t.corr);     // Antisymmetric filter
    }
  }
};

#endif

