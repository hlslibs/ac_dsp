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

//*****************************************************************************************//
//
// File: ac_fir_prog_coeffs.h
//
// Description:
// This file contains wrapper classes that implement different architectures of FIR
// filters. Each top level class calls a different member function of the core
// class and hence implements a different architecture.
// The fir_prog_coeffs_core is core class of the filter and has filter functions
// definitions.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_fir_prog_coeffs.h>
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Define input, output, coefficient and accumulator type
//      typedef ac_fixed<32, 16, true> IN_TYPE;
//      typedef ac_fixed<64, 32, true> OUT_TYPE;
//      typedef ac_fixed<32, 16, true> COEFF_TYPE;
//      typedef ac_fixed<64, 32, true> MAC_TYPE;
//      ac_fir_load_coeffs< IN_TYPE, OUT_TYPE, COEFF_TYPE, MAC_TYPE, N_TAPS, FOLD_ODD > filter_design1;
//      // Initialize channels for input, output, coefficients and load flag
//      ac_channel<IN_TYPE>  input;
//      ac_channel<OUT_TYPE> output;
//      ac_channel<COEFF_TYPE> coeffs_ch;
//      // Make sure to initialize the inputs and the coeffs
//      // Call the top-level function.
//      filter_design1.run(input, coeffs_ch, output, load);
//
//      CCS_RETURN(0);
//    }
//
//***************************************************************************************//

#ifndef _INCLUDED_AC_FIR_PROG_COEFFS_H_
#define _INCLUDED_AC_FIR_PROG_COEFFS_H_

#include <ac_fixed.h>
#include <ac_int.h>

// Make sure that the enum is only defined once.
#ifndef __FIR_FILTER_TYPES_ENUM_DEF__
#define __FIR_FILTER_TYPES_ENUM_DEF__
// The parameters within this enum help the user choose between different filter architectures.
typedef enum { SHIFT_REG, ROTATE_SHIFT, C_BUFF, FOLD_EVEN, FOLD_ODD, TRANSPOSED, FOLD_EVEN_ANTI, FOLD_ODD_ANTI } FTYPE;

#endif



/*
 * class Declarations
 */

template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, int N_TAPS >
class fir_prog_coeffs_core
{
private: // Data
  IN_TYPE reg[N_TAPS];                                      // shift register for filter implementation
  ACC_TYPE reg_trans[N_TAPS];                               // register for transpose form filter implementation
  ac_int < ac::log2_ceil < N_TAPS >::val + 1, false > wptr; // Read pointer for circular buffer implementation
  ac_int < ac::log2_ceil < N_TAPS >::val + 1, true > rptr;  // write pointer for circular buffer implementation

public: // Functions
  // Constructor
  fir_prog_coeffs_core() {
    ac::init_array < AC_VAL_0 > (reg, N_TAPS);
    ac::init_array < AC_VAL_0 > (reg_trans, N_TAPS);
    wptr = 0;
    rptr = 0;
  }

  // firShiftReg() implements shift register for the Filter
  void firShiftReg(IN_TYPE din) {
#pragma hls_unroll yes
    SHIFT:
    for (int i = (N_TAPS - 1); i >= 0; i--) {
      reg[i] = (i == 0) ? din : reg[i - 1];
    }
  }

  // firCircularBuffWrite() implements a circular buffer write operation
  void firCircularBuffWrite(IN_TYPE din) {
    reg[wptr] = din;
    if (wptr == N_TAPS - 1) {
      wptr = 0;
    } else {
      wptr++;
    }
  }

  // firCircularBuffRead() implements a circular buffer read operation
  IN_TYPE firCircularBuffRead(ac_int < ac::log2_ceil < N_TAPS >::val, false > idx) {
    rptr = wptr - 1 - idx;
    if (rptr < 0) {
      rptr = rptr + N_TAPS;
    }
    return reg[rptr];
  }

  // firProgCoeffsShiftReg() Implements a non symmetric FIR filter with shift register
  void firProgCoeffsShiftReg(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    firShiftReg(data_in);
    MAC:
    for (int i = (N_TAPS - 1); i >= 0; i--) {
      acc += reg[i] * coeffs[i];
    }
    data_out = acc;
  }

  // firProgCoeffsRotateShift() implements a rotational shift based implementation
  void firProgCoeffsRotateShift(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    IN_TYPE temp_rotate;
    MAC:
    for (int i = (N_TAPS); i >= 0; i--) {
      if (i == N_TAPS) {
        temp_rotate = data_in;
      } else {
        temp_rotate = reg[N_TAPS - 1];
        acc += reg[N_TAPS - 1] * coeffs[i];
      }
      firShiftReg(temp_rotate);
    }
    data_out = acc;
  }

  // firProgCoeffsCircularBuff() Circular Buffer based implementation
  void firProgCoeffsCircularBuff(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    MAC:
    for (int i = 0; i <= (N_TAPS - 1); i++) {
      if (i == 0) {
        firCircularBuffWrite(data_in);      // Store in circular buffer
      }
      acc += firCircularBuffRead(i) * coeffs[i];   // Read from circular buffer
    }
    data_out = acc;
  }

  // firProgCoeffsShiftRegSymmetricEvenTaps() implements a symmetric filter with even number of Taps
  // and shift register based implementation
  void firProgCoeffsShiftRegSymmetricEvenTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    firShiftReg(data_in);
    MAC:
    for (int i = (N_TAPS / 2) - 1; i >= 0; i--) {
      acc += coeffs[i] * (reg[i] + reg[N_TAPS - 1 - i]);
    }
    data_out = acc;
  }

  // firProgCoeffsShiftRegSymmetricOddTaps() implements a symmetric filter with odd number of Taps
  // and shift register based implementation
  void firProgCoeffsShiftRegSymmetricOddTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    ACC_TYPE fold = 0;
    firShiftReg(data_in);
    MAC:
    for (int i = 0; i < (((N_TAPS - 1) / 2) + 1); i++) {
      if (i == (N_TAPS - 1) / 2) {
        fold = reg[i]; // Central Tap passes through
      } else {
        fold = reg[i] + reg[(N_TAPS - 1) - i];
      }
      acc += coeffs[i] * fold;
    }
    data_out = acc;
  }

  // firProgCoeffsTransposed() implements a transpose form of the filter
  void firProgCoeffsTransposed(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE temp = 0;
    IN_TYPE in = data_in;
#pragma unroll yes
    MAC:
    for (int i = (N_TAPS - 1); i >= 0; i--) {
      if (i == 0) {
        temp = 0;
      } else {
        temp = reg_trans[i - 1];
      }
      reg_trans[i] = in * coeffs[(N_TAPS - 1) - i] + temp;
    }
    data_out = reg_trans[N_TAPS - 1];
  }

};

#include <ac_channel.h>

template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, int N_TAPS, FTYPE ftype = SHIFT_REG >
class ac_fir_prog_coeffs
{
private: // Data
  fir_prog_coeffs_core < IN_TYPE, OUT_TYPE, COEFF_TYPE, ACC_TYPE, N_TAPS > filter;

  // Internal array that stores coefficients data
  COEFF_TYPE coeffs[N_TAPS];

public: // Functions
  // Constructor
  ac_fir_prog_coeffs() {
    ac::init_array < AC_VAL_DC > (coeffs, N_TAPS);
  };
#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  // run() is top function for C++ module. Based on filter type configured
  // it instantiates the desired filter function
  void run(ac_channel < IN_TYPE > &data_in, ac_channel < COEFF_TYPE > &coeffs_ch, ac_channel < OUT_TYPE > &data_out) {
    IN_TYPE core_in;
    OUT_TYPE core_out;

    if (coeffs_ch.available(N_TAPS)) {
      for (int i = 0; i < N_TAPS; i++) { coeffs[i] = coeffs_ch.read(); }
    }

    if (data_in.available(1)) {
      core_in = data_in.read();
      if (ftype == SHIFT_REG) {
        filter.firProgCoeffsShiftReg(core_in, coeffs, core_out);                  // Shift register implementation
      }
      if (ftype == ROTATE_SHIFT) {
        filter.firProgCoeffsRotateShift(core_in, coeffs, core_out);               // Rotate shift implementation
      }
      if (ftype == C_BUFF) {
        filter.firProgCoeffsCircularBuff(core_in, coeffs, core_out);              // Circular buffer based implementation
      }
      if (ftype == FOLD_EVEN) {
        filter.firProgCoeffsShiftRegSymmetricEvenTaps(core_in, coeffs, core_out); // Symmetric filter with even number of Taps
      }
      if (ftype == FOLD_ODD) {
        filter.firProgCoeffsShiftRegSymmetricOddTaps(core_in, coeffs, core_out);  // Symmetric filter with odd number of Taps
      }
      if (ftype == TRANSPOSED) {
        filter.firProgCoeffsTransposed(core_in, coeffs, core_out);                // Transposed form of the filter
      }
      data_out.write(core_out);
    }
  }

};

#endif

