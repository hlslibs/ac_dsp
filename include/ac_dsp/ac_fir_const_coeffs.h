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
// File: ac_fir_const_coeffs.h
//
// Description:
//    This file contains the C++ interface to the core FIR class, as well as the core FIR class which
//    contains the core functionality of the filter. The member function "run()" contained in the
//    ac_fir_const_coeffs class is the top level function for the C++ design.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_fir_const_coeffs.h>
//    #include <mc_scverify.h>
//
//    #pragma hls_design top
//    template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, unsigned N_TAPS, FTYPE ftype >
//    class ac_fir_const_coeffs_wrapper : public ac_fir_const_coeffs<IN_TYPE, OUT_TYPE, COEFF_TYPE, ACC_TYPE, N_TAPS, ftype>
//    {
//    private:
//      const COEFF_TYPE coeffs[N_TAPS] = {1, 2, 3, 4, 5};
//    public:
//      ac_fir_const_coeffs_wrapper() : ac_fir_const_coeffs<IN_TYPE, OUT_TYPE, COEFF_TYPE, ACC_TYPE, N_TAPS, ftype> (coeffs) {  }
//    };
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
//      // Make sure that this number of taps is consistent with the wrapper class coeffs array.
//      ac_fir_const_coeffs_wrapper< IN_TYPE_TB, OUT_TYPE_TB, COEFF_TYPE_TB, MAC_TYPE_TB, 5, FOLD_ODD > filter_design1;
//      // Initialize channels for input and output
//      ac_channel<IN_TYPE_TB>  input;
//      ac_channel<OUT_TYPE_TB> output;
//      // Write 4 inputs to the input port
//      input.write(.125);
//      input.write(.375);
//      input.write(1.525);
//      input.write(.125);
//      // Call the top-level function.
//      filter_design1.run(input, output);
//
//      CCS_RETURN(0);
//    }
//
//*********************************************************************************************************

#ifndef _INCLUDED_AC_FIR_CONST_COEFFS_H_
#define _INCLUDED_AC_FIR_CONST_COEFFS_H_

#include <ac_fixed.h>
#include <ac_int.h>
#include <ac_channel.h>

//****************************************************************************************
//
// Class : fir_const_coeffs_core
//
// Description:
//   class "fir_const_coeffs_core" is the core class for constant coefficient FIR filters.
//   The class member functions implement different architectures for FIR filter.
//
//****************************************************************************************


template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, unsigned N_TAPS >
class fir_const_coeffs_core
{
private: // Data
  const COEFF_TYPE *const coeffs;

  IN_TYPE reg[N_TAPS];                                       // shift register for filter implementation
  ACC_TYPE reg_trans[N_TAPS];                                // register for transpose form filter implementation
  ac_int < ac::log2_ceil < N_TAPS >::val + 1, false > wptr;  // Read pointer for circular buffer implementation
  ac_int < ac::log2_ceil < N_TAPS >::val + 1, true > rptr;   // write pointer for circular buffer implementation

public: // Functions
  // Constructor
  fir_const_coeffs_core(const COEFF_TYPE *const ptr) : coeffs(ptr) {
    wptr = 0;
    rptr = 0;
    ac::init_array < AC_VAL_0 > (reg, N_TAPS);
    ac::init_array < AC_VAL_0 > (reg_trans, N_TAPS);
  }

  // firShiftReg() implements a shift register for the Filter
  void firShiftReg(IN_TYPE din) {
#pragma hls_unroll yes
    SHIFT:
    for (int i = (N_TAPS - 1); i >= 0; i--) {
      reg[i] = (i == 0) ? din : reg[i - 1];
    }
  }

  // firCircularBuffWrite() implements a Circular buffer write operation
  void firCircularBuffWrite(IN_TYPE din) {
    reg[wptr] = din;
    if (wptr == N_TAPS - 1) {
      wptr = 0;
    } else {
      wptr++;
    }
  }

  // firCircularBuffRead() implements a Circular buffer read operation
  IN_TYPE firCircularBuffRead(ac_int < ac::log2_ceil < N_TAPS >::val, false > idx) {
    rptr = wptr - 1 - idx;
    if (rptr < 0) {
      rptr = rptr + N_TAPS;
    }
    return reg[rptr];
  }

  // firConstCoeffsShiftReg() implements a non symmetric FIR filter with shift register
  void firConstCoeffsShiftReg(IN_TYPE &data_in, OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    firShiftReg(data_in);       // Shift Reg
    MAC:
    for (int i = (N_TAPS - 1); i >= 0; i--) {
      acc += reg[i] * coeffs[i];
    }
    data_out = acc;
  }

  // firConstCoeffsRotateShift() implements a Rotational shift based implementation
  void firConstCoeffsRotateShift(IN_TYPE &data_in, OUT_TYPE &data_out) {
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
      firShiftReg(temp_rotate);       // Rotate right
    }
    data_out = acc;
  }

  // firConstCoeffsCircularBuff() implements a Circular Buffer based implementation
  void firConstCoeffsCircularBuff(IN_TYPE &data_in, OUT_TYPE &data_out) {
    ACC_TYPE acc = 0.0;
    MAC:
    for (int i = 0; i <= (N_TAPS - 1); i++) {
      if (i == 0) {
        firCircularBuffWrite(data_in);              // Store in circular buffer
      }
      acc += firCircularBuffRead(i) * coeffs[i];    // Read from circular buffer
    }
    data_out = acc;
  }

  // firConstCoeffsShiftRegSymmetricEvenTaps() implements a symmetric filter with even number of Taps
  // and shift register based implementation
  void firConstCoeffsShiftRegSymmetricEvenTaps(IN_TYPE &data_in, OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    firShiftReg(data_in);
    MAC:
    for (int i = (N_TAPS / 2) - 1; i >= 0; i--) {
      acc += coeffs[i] * (reg[i] + reg[N_TAPS - 1 - i]);
    }
    data_out = acc;
  }

  // firConstCoeffsShiftRegSymmetricOddTaps() implements a symmetric filter with odd number of Taps
  // and shift register based implementation
  void firConstCoeffsShiftRegSymmetricOddTaps(IN_TYPE &data_in, OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    ACC_TYPE fold = 0;
    firShiftReg(data_in);
    MAC:
    for (int i = 0; i < (((N_TAPS - 1) / 2) + 1); i++) {
      if (i == (N_TAPS - 1) / 2) {
        fold = reg[i];  // Central tap passes through
      } else {
        fold = reg[i] + reg[(N_TAPS - 1) - i];
      }
      acc += coeffs[i] * fold;
    }
    data_out = acc;
  }

  // firConstCoeffsTransposed() implements the transposed form of the filter
  void firConstCoeffsTransposed(IN_TYPE &data_in, OUT_TYPE &data_out) {
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

// Make sure that this enum is only defined once, as different FIR designs use and define the same enum.
#ifndef __FIR_FILTER_TYPES_ENUM_DEF__
#define __FIR_FILTER_TYPES_ENUM_DEF__
// The parameters within this enum help the user choose between different filter architectures.
typedef enum { SHIFT_REG, ROTATE_SHIFT, C_BUFF, FOLD_EVEN, FOLD_ODD, TRANSPOSED, FOLD_EVEN_ANTI, FOLD_ODD_ANTI } FTYPE;
#endif

//********************************************************************************************************************
//
// Class : ac_fir_const_coeffs
// Description:
//   This class contains the top function for the constant coefficient FIR filter ('run()'). It has the pointer to the
//   coeffs array passed as a constructor parameter from the wrapper class, which derives it.
//
//********************************************************************************************************************

template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, unsigned N_TAPS, FTYPE ftype >
class ac_fir_const_coeffs
{
private:
  fir_const_coeffs_core < IN_TYPE, OUT_TYPE, COEFF_TYPE, ACC_TYPE, N_TAPS > filter;

public:

  // constructor with pointer of const coeff array as an arg
  ac_fir_const_coeffs(const COEFF_TYPE *const c_ptr) : filter(c_ptr) { }

#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void run(ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out) {
    IN_TYPE core_in;
    OUT_TYPE core_out;
#ifndef __SYNTHESIS__
    while (data_in.available(1))
#endif
    {
      core_in = data_in.read();
      if (ftype == SHIFT_REG) {
        filter.firConstCoeffsShiftReg(core_in, core_out);                   // Shift register implementation
      }
      if (ftype == ROTATE_SHIFT) {
        filter.firConstCoeffsRotateShift(core_in, core_out);                // Rotate shift implementation
      }
      if (ftype == C_BUFF) {
        filter.firConstCoeffsCircularBuff(core_in, core_out);               // Circular buffer based implementation
      }
      if (ftype == FOLD_EVEN) {
        filter.firConstCoeffsShiftRegSymmetricEvenTaps(core_in, core_out);  // Symmetric filter with even number of Taps
      }
      if (ftype == FOLD_ODD) {
        filter.firConstCoeffsShiftRegSymmetricOddTaps(core_in, core_out);   // Symmetric filter with odd number of Taps
      }
      if (ftype == TRANSPOSED) {
        filter.firConstCoeffsTransposed(core_in, core_out);                 // Transposed form of the filter
      }
      data_out.write(core_out);
    }
  }
};

#endif

