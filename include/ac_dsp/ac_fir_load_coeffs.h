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

//*********************************************************************************************************
//
// File: ac_fir_load_coeffs.h
//
// Description:
//    This file contains the C++ interface to the core FIR class, as well as the core FIR class which
//    contains the core functionality of the filter. The member function "run()" contained in the
//    ac_fir_load_coeffs class is the top level function for the C++ design. The core class is generic,
//    and can be used by both the supplied C++ implementation as well as any other SystemC implementations
//    the user might want to use.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_fir_load_coeffs.h>
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Define input, output, coefficient and accumulator type
//      typedef ac_fixed<32, 16, true> IN_TYPE_TB;
//      typedef ac_fixed<64, 32, true> OUT_TYPE_TB;
//      typedef ac_fixed<32, 16, true> COEFF_TYPE_TB;
//      typedef ac_fixed<64, 32, true> MAC_TYPE_TB;
//      // We are using the folded version of a symmetric filter with an odd number of taps for this design.
//      // Hence, we pass FOLD_ODD as an additional template parameter. The filter itself has 5 taps.
//      ac_fir_load_coeffs< IN_TYPE_TB, OUT_TYPE_TB, COEFF_TYPE_TB, MAC_TYPE_TB, 5, FOLD_ODD > filter_design1;
//      // Initialize channels for input, output, coefficients and load flag
//      ac_channel<IN_TYPE_TB>  input;
//      ac_channel<OUT_TYPE_TB> output;
//      ac_channel<COEFF_TYPE_TB> coeffs_ch;
//      ac_channel<bool> load;
//      // Write 5 coefficients to the coefficient channel
//      coeffs_ch.write(2.5);
//      coeffs_ch.write(1.25);
//      coeffs_ch.write(0.5);
//      coeffs_ch.write(0.375);
//      coeffs_ch.write(1.75);
//      // Write the value "true" to the coeffs_ch channel.
//      load.write(true);
//      // Call the top-level function to load all the coeff values.
//      filter_design1.run(input, coeffs_ch, output, load);
//      // Write the value "false" to the coeffs_ch channel to signal the design to stop loading coeffs values.
//      load.write(false);
//      // Write 5 inputs to the input port
//      input.write(.125);
//      input.write(.375);
//      input.write(1.525);
//      input.write(.125);
//      input.write(.75);
//      // Call the top-level function again to process the input values.
//      filter_design1.run(input, coeffs_ch, output, load);
//
//      CCS_RETURN(0);
//    }
//
//*********************************************************************************************************

#ifndef _INCLUDED_AC_FIR_LOAD_COEFFS_H_
#define _INCLUDED_AC_FIR_LOAD_COEFFS_H_

#include <ac_fixed.h>
#include <ac_int.h>

#ifndef __SYNTHESIS__
#include <iostream>
using namespace std;
#endif

//****************************************************************************************
//
// Class : fir_load_coeffs_core
//
// Description:
//   class "fir_load_coeffs_core" is the core class for loadable coefficient FIR filters.
//   The class member functions implement different architectures for FIR filter.
//
//****************************************************************************************

template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, unsigned N_TAPS >
class fir_load_coeffs_core
{

private: // Data
  IN_TYPE reg[N_TAPS];                                      // shift register for filter implementation
  ACC_TYPE reg_trans[N_TAPS];                               // register for transpose form filter implementation
  ac_int < ac::log2_ceil < N_TAPS >::val + 1, false > wptr; // Read pointer for circular buffer implementation
  ac_int < ac::log2_ceil < N_TAPS >::val + 1, true > rptr;  // write pointer for circular buffer implementation

public: // Functions
  // Constructor
  fir_load_coeffs_core() {
    ac::init_array < AC_VAL_0 > (reg, N_TAPS);
    ac::init_array < AC_VAL_0 > (reg_trans, N_TAPS);
    wptr = 0;
    rptr = 0;
  }

  // firShiftReg() implements a shift register for the Filter
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
    if (wptr == N_TAPS - 1)
    { wptr = 0; }
    else
    { wptr++; }
  }

  // firCircularBuffRead() implements a circular buffer read operation
  IN_TYPE firCircularBuffRead(ac_int < ac::log2_ceil < N_TAPS >::val, false > idx) {
    rptr = wptr - 1 - idx;
    if (rptr < 0)
    { rptr = rptr + N_TAPS; }
    return reg[rptr];
  }

  // firLoadCoeffsShiftReg() Implements a non symmetric FIR filter with shift register
  void firLoadCoeffsShiftReg(IN_TYPE &data_in, COEFF_TYPE h[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    firShiftReg(data_in);
    MAC:
    for (int i = (N_TAPS - 1); i >= 0; i--) {
      acc += reg[i] * h[i];
    }
    data_out = acc;
  }

  // firLoadCoeffsRotateShift() implements a rotational shift based implementation
  void firLoadCoeffsRotateShift(IN_TYPE &data_in, COEFF_TYPE h[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    IN_TYPE temp_rotate;
    MAC:
    for (int i = (N_TAPS); i >= 0; i--) {
      if (i == N_TAPS) {
        temp_rotate = data_in;
      } else {
        temp_rotate = reg[N_TAPS - 1];
        acc += reg[N_TAPS - 1] * h[i];
      }
      firShiftReg(temp_rotate);
    }
    data_out = acc;
  }

  // firLoadCoeffsCircularBuff() implements a Circular Buffer based implementation
  void firLoadCoeffsCircularBuff(IN_TYPE &data_in, COEFF_TYPE h[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    MAC:
    for (int i = 0; i <= (N_TAPS - 1); i++) {
      if (i == 0) {
        firCircularBuffWrite(data_in);      // Store in circular buffer
      }
      acc += firCircularBuffRead(i) * h[i];   // Read from circular buffer
    }
    data_out = acc;
  }

  // firLoadCoeffsShiftRegSymmetricEvenTaps() implements a symmetric filter with even number of Taps
  // and shift register based implementation
  void firLoadCoeffsShiftRegSymmetricEvenTaps(IN_TYPE &data_in, COEFF_TYPE h[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    firShiftReg(data_in);
    MAC:
    for (int i = (N_TAPS / 2) - 1; i >= 0; i--) {
      acc += h[i] * (reg[i] + reg[N_TAPS - 1 - i]);
    }
    data_out = acc;
  }

  // firLoadCoeffsShiftRegSymmetricOddTaps() implements a symmetric filter with odd number of Taps
  // and shift register based implementation
  void firLoadCoeffsShiftRegSymmetricOddTaps(IN_TYPE &data_in, COEFF_TYPE h[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    ACC_TYPE fold = 0;
    firShiftReg(data_in);
    MAC:
    for (int i = 0; i < (((N_TAPS - 1) / 2) + 1); i++) {
      if (i == (N_TAPS - 1) / 2)
      { fold = reg[i]; }      // Central Tap passes through
      else
      { fold = reg[i] + reg[(N_TAPS - 1) - i]; }
      acc += h[i] * fold;
    }
    data_out = acc;
  }

  // firLoadCoeffsTransposed() implements a transposed form of the filter
  void firLoadCoeffsTransposed(IN_TYPE &data_in, COEFF_TYPE h[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE temp = 0;
    IN_TYPE in = data_in;
#pragma unroll yes
    MAC:
    for (int i = (N_TAPS - 1); i >= 0; i--) {
      if (i == 0)
      { temp = 0; }
      else
      { temp = reg_trans[i - 1]; }
      reg_trans[i] = in * h[(N_TAPS - 1) - i] + temp;
    }
    data_out = reg_trans[N_TAPS - 1];
  }

};

#include <ac_channel.h>

// Make sure that this enum is only defined once, as different FIR designs use and define the same enum.
#ifndef __FIR_FILTER_TYPES_ENUM_DEF__
#define __FIR_FILTER_TYPES_ENUM_DEF__
// The parameters within this enum help the user choose between different filter architectures.
typedef enum { SHIFT_REG, ROTATE_SHIFT, C_BUFF, FOLD_EVEN, FOLD_ODD, TRANSPOSED, FOLD_EVEN_ANTI, FOLD_ODD_ANTI } FTYPE;
#endif

//*******************************************************************************************
//
// Class : ac_fir_load_coeffs
// Description:
//   This class contains the top function for the loadable coefficients FIR filter ('run()').
//
//*******************************************************************************************

// Top level class for the Shift register implementation
template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, unsigned N_TAPS, FTYPE ftype >
class ac_fir_load_coeffs
{
private: // Data
  fir_load_coeffs_core < IN_TYPE, OUT_TYPE, COEFF_TYPE, ACC_TYPE, N_TAPS > filter;
  // Internal array that stores coefficient data
  COEFF_TYPE coeffs[N_TAPS];

public: // Functions
  // Constructor
  ac_fir_load_coeffs() {
    ac::init_array < AC_VAL_DC > (coeffs, N_TAPS);
  };

#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void run(ac_channel < IN_TYPE > &data_in, ac_channel <COEFF_TYPE> &coeffs_ch, ac_channel < OUT_TYPE > &data_out, ac_channel< bool > &ld) {
    IN_TYPE core_in;
    OUT_TYPE core_out;

    if (ld.available(1)) {
      bool ld_t;
      ld_t = ld.read();
      // Whenever the load flag is set to true and all coefficients are available on the coeffs_ch channel, load the coefficient values into the coeffs array
      if (ld_t && coeffs_ch.available(N_TAPS)) {
        for (int i = 0; i < N_TAPS; i++) { coeffs[i] = coeffs_ch.read(); }
      }
    }


#ifndef __SYNTHESIS__
    while (data_in.available(1))
#endif
    {
      core_in = data_in.read();
      if (ftype == SHIFT_REG) {
        filter.firLoadCoeffsShiftReg(core_in, coeffs, core_out); // Shift register implementation
      }
      if (ftype == ROTATE_SHIFT) {
        filter.firLoadCoeffsRotateShift(core_in, coeffs, core_out); // Rotate shift implementation
      }
      if (ftype == C_BUFF) {
        filter.firLoadCoeffsCircularBuff(core_in, coeffs, core_out); // Circular buffer based implementation
      }
      if (ftype == FOLD_EVEN) {
        filter.firLoadCoeffsShiftRegSymmetricEvenTaps(core_in, coeffs, core_out); // Symmetric filter with even number of Taps
      }
      if (ftype == FOLD_ODD) {
        filter.firLoadCoeffsShiftRegSymmetricOddTaps(core_in, coeffs, core_out); // Symmetric filter with odd number of Taps
      }
      if (ftype == TRANSPOSED) {
        filter.firLoadCoeffsTransposed(core_in, coeffs, core_out); // Transposed form of the filter
      }
      data_out.write(core_out);
    }
  }
};

#endif

