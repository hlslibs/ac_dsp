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

//************************************************************************************
// File:         ac_fir_reg_share.h
//
// Description:
// 		This file contains wrapper classes that implement different architectures of FIR
// 		filters. Each top level class calls a different member function of the core
// 		class and hence implements a different architecture.
// 		The fir_reg_share_core is core class of the filter and has filter functions
// 		definitions
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
//      ac_fir_reg_share< IN_TYPE, OUT_TYPE, COEFF_TYPE, MAC_TYPE, N_TAPS, MEM_WORD_WIDTH,
//      BLK_SZ, BLK_OFFSET, FOLD_ODD > filter_design1;
//      IN_TYPE  input;
//      OUT_TYPE output;
//      COEFF_TYPE coeffs;
//      // Make sure to initialize the inputs and the coeffs
//      // Call the top-level function.
//      filter_design1.run(input, coeffs, output);
//
//      CCS_RETURN(0);
//    }
//
//*************************************************************************************

#ifndef _INCLUDED_AC_FIR_REG_SHARE_H_
#define _INCLUDED_AC_FIR_REG_SHARE_H_

#include <ac_fixed.h>
#include <ac_int.h>

// Make sure that the enum is only defined once.
#ifndef __FIR_FILTER_TYPES_ENUM_DEF__
#define __FIR_FILTER_TYPES_ENUM_DEF__
// The parameters within this enum help the user choose between different filter architectures.
typedef enum { SHIFT_REG, ROTATE_SHIFT, C_BUFF, FOLD_EVEN, FOLD_ODD, TRANSPOSED, FOLD_EVEN_ANTI, FOLD_ODD_ANTI } FTYPE;

#endif

typedef ac_fixed < 16,1, true > DEFAULT_TYPE;

//class Declarations


template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, int N_TAPS, int MEM_WORD_WIDTH, int BLK_SZ, int BLK_OFFSET >
class fir_reg_share_core
{

private: // Data
  IN_TYPE *reg;

public: // Functions
  // Constructor
  fir_reg_share_core() { };

  void fir_reg_share_core_init(IN_TYPE *ptr_t) {
    reg = ptr_t;
  }

  // firShiftReg() implements shift register for the Filter
  void firShiftReg(IN_TYPE din) {
#pragma hls_unroll yes
    SHIFT: for (int i = (N_TAPS - 1); i >= 0; i--) {
      reg[i] = (i == 0) ? din : reg[i - 1];
    }
  }

  void firProgCoeffsDelayLine(OUT_TYPE &data_out) {
    data_out = reg[N_TAPS -1];
  }

  // firProgCoeffsShiftReg() Implements non symmetric FIR filter with shift register
  void firProgCoeffsShiftReg(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    int ram_addr = 0;
    int index = 0;
    BLK: for ( int i = 0; i < N_TAPS; i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH )  {
      index = 0;
#pragma unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
        acc += (reg[i + index]) * coeffs[ram_addr + block_count];
        index++;
      }
    }
    data_out = acc;
  }

  // firProgCoeffsShiftRegSymmetricEvenTaps() symmetric filter with even number of Taps
  // and shift register based implementation
  void firProgCoeffsShiftRegSymmetricEvenTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    int ram_addr = 0;
    int index = 0;
    BLK: for ( int i = 0; i < N_TAPS/2; i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH )  {
      index = 0;
#pragma unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
        acc += (reg[i + index] + reg[N_TAPS - 1 - i - index] ) * coeffs[ram_addr + block_count];
        index++;
      }
    }
    data_out = acc;
  }

  void firProgCoeffsShiftRegAntiSymmetricEvenTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    int ram_addr = 0;
    int index = 0;
    BLK: for ( int i = 0; i < N_TAPS/2; i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH )  {
      index = 0;
#pragma unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
        acc += (reg[i + index] - reg[N_TAPS - 1 - i - index] ) * coeffs[ram_addr + block_count];
        index++;
      }
    }
    data_out = acc;
  }

  // firProgCoeffsShiftRegSymmetricOddTaps() symmetric filter with odd number of Taps
  // and shift register based implementation
  void firProgCoeffsShiftRegSymmetricOddTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    ACC_TYPE fold = 0;
    int ram_addr = 0;
    int index = 0;
    BLK: for (int i = 0; i < (((N_TAPS - 1) / 2) + 1); i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH) {
      index = 0;
#pragma unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
        if ((i + index)  == (N_TAPS - 1) / 2)
        { fold = reg[i + index]; }                      // Central passes through
        else
        { fold = reg[i + index] + reg[(N_TAPS - 1) - i - index]; }
        index++;
        acc += coeffs[ram_addr + block_count] * fold;
      }
    }
    data_out = acc;
  }

  void firProgCoeffsShiftRegAntiSymmetricOddTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0;
    ACC_TYPE fold = 0;
    int ram_addr = 0;
    int index = 0;
    BLK: for (int i = 0; i < (((N_TAPS - 1) / 2) + 1); i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH) {
      index = 0;
#pragma unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
        if ((i + index)  == (N_TAPS - 1) / 2)
        { fold = reg[i + index]; }                      // Central passes through
        else
        { fold = reg[i + index] - reg[(N_TAPS - 1) - i - index]; }
        index++;
        acc += coeffs[ram_addr + block_count] * fold;
      }
    }
    data_out = acc;
  }

};

template < int N_TAPS = 2, class IN_TYPE = DEFAULT_TYPE, class OUT_TYPE = DEFAULT_TYPE, class COEFF_TYPE = DEFAULT_TYPE, class ACC_TYPE = DEFAULT_TYPE, int MEM_WORD_WIDTH = 1, int BLK_SZ = 1, int BLK_OFFSET = 0,FTYPE ftype = SHIFT_REG >
class ac_fir_reg_share
{
private: // Data
  IN_TYPE *ptr;
  fir_reg_share_core < IN_TYPE, OUT_TYPE, COEFF_TYPE, ACC_TYPE, N_TAPS, MEM_WORD_WIDTH, BLK_SZ, BLK_OFFSET > filter;
  // Internal array that stores coefficients data
  COEFF_TYPE coeffs[N_TAPS];

public: // Functions
  // Constructor
  ac_fir_reg_share(IN_TYPE *ptr_t) {
    ptr = ptr_t;
  };

  // Based on filter type configured it instantiates the desired filter function
  void run(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    IN_TYPE core_in;
    OUT_TYPE core_out;

    filter.fir_reg_share_core_init(ptr);
    filter.firShiftReg(data_in);
    core_in = data_in;
    filter.fir_reg_share_core_init(ptr);
    if (ftype == SHIFT_REG) {
      filter.firProgCoeffsShiftReg(core_in, coeffs, core_out); // Shift register implementation
    }
    if (ftype == FOLD_EVEN) {
      filter.firProgCoeffsShiftRegSymmetricEvenTaps(core_in, coeffs, core_out); // Symmetric filter with even number of Taps
    }
    if (ftype == FOLD_EVEN_ANTI) {
      filter.firProgCoeffsShiftRegAntiSymmetricEvenTaps(core_in, coeffs, core_out); // Anti Symmetric filter with even number of Taps
    }
    if (ftype == FOLD_ODD) {
      filter.firProgCoeffsShiftRegSymmetricOddTaps(core_in, coeffs, core_out); // Symmetric filter with odd number of Taps
    }
    if (ftype == FOLD_ODD_ANTI) {
      filter.firProgCoeffsShiftRegAntiSymmetricOddTaps(core_in, coeffs, core_out); // Anti Symmetric filter with odd number of Taps
    }
    data_out = core_out;
  }

  void ac_firProgCoeffs_delay_line(OUT_TYPE &core_out) {
    filter.fir_reg_share_core_init(ptr);
    filter.firProgCoeffsDelayLine(core_out);
  }

};

#endif

