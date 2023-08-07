/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) DSP Library                                        *
 *                                                                        *
 *  Software Version: 3.5                                                 *
 *                                                                        *
 *  Release Date    : Sun Jul 23 16:34:46 PDT 2023                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.5.0                                               *
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
// Revision History:
//    3.4.0  - Added CDesignChecker waivers/fixes for ac_dsp IP blocks.
//             Changes made in general:
//               - CNS violations were waived away.
//               - FXD violations were fixed by typecasting and assigning fixed value
//                 initializations to ac_fixed variables.
//               - ABR and ABW violations were waived away.
//
//
//*************************************************************************************

#ifndef _INCLUDED_AC_FIR_REG_SHARE_H_
#define _INCLUDED_AC_FIR_REG_SHARE_H_

#include <ac_fixed.h>
#include <ac_int.h>
#include <ac_channel.h>

// Make sure that the enum is only defined once.
#ifndef __FIR_FILTER_TYPES_ENUM_DEF__
#define __FIR_FILTER_TYPES_ENUM_DEF__
// The parameters within this enum help the user choose between different filter architectures.
typedef enum { SHIFT_REG, ROTATE_SHIFT, C_BUFF, FOLD_EVEN, FOLD_ODD, TRANSPOSED, FOLD_EVEN_ANTI, FOLD_ODD_ANTI } FTYPE;

#endif

#include <mc_scverify.h>

typedef ac_fixed < 16,1, true > DEFAULT_TYPE;

//class Declarations

//===================================================================================================================
// Class: fir_reg_share_core
// Description: The class member functions implement different architectures for FIR filter.
//-------------------------------------------------------------------------------------------------------------------

template < class IN_TYPE, class OUT_TYPE, class COEFF_TYPE, class ACC_TYPE, int N_TAPS, int MEM_WORD_WIDTH, int BLK_SZ, int BLK_OFFSET >
class fir_reg_share_core
{

private: // Data
  IN_TYPE *reg;

public: 
  // Constructor
  fir_reg_share_core() { };

  void fir_reg_share_core_init(IN_TYPE *ptr_t) {
    reg = ptr_t;
  }

//------------------------------------------------------------------------------------------------------
// Member Function: firShiftReg()
// Description: firShiftReg() implements a shift register for the Filter

  void firShiftReg(IN_TYPE din) {
#pragma hls_unroll yes
    SHIFT: for (int i = (N_TAPS - 1); i >= 0; i--) {
#pragma hls_waive ABR ABW
      reg[i] = (i == 0) ? din : reg[i - 1];
    }
  }

//------------------------------------------------------------------------------------------------------
// Member Function: firProgCoeffsDelayLine()
// Description: firProgCoeffsDelayLine() implements a delay line

  void firProgCoeffsDelayLine(OUT_TYPE &data_out) {
    data_out = reg[N_TAPS -1];
  }

//------------------------------------------------------------------------------------------------------
// Member Function: firProgCoeffsShiftReg()
// Description: firProgCoeffsShiftReg() implements a non symmetric FIR filter with shift register

  void firProgCoeffsShiftReg(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0.0;
    int ram_addr = 0;
    int index = 0;
    BLK: for ( int i = 0; i < N_TAPS; i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH )  {
      index = 0;
#pragma hls_unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
#pragma hls_waive ABR
        acc += (reg[i + index]) * coeffs[ram_addr + block_count];
        index++;
      }
    }
    data_out = acc;
  }

//----------------------------------------------------------------------------------------------------------------
// Member Function: firProgCoeffsShiftRegSymmetricEvenTaps()
// Description: firProgCoeffsShiftRegSymmetricEvenTaps() implements a symmetric filter with even number of Taps
// and shift register based implementation

  void firProgCoeffsShiftRegSymmetricEvenTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0.0;
    int ram_addr = 0;
    int index = 0;
    BLK: for ( int i = 0; i < N_TAPS/2; i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH )  {
      index = 0;
#pragma hls_unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
#pragma hls_waive ABR
        acc += (reg[i + index] + reg[N_TAPS - 1 - i - index] ) * coeffs[ram_addr + block_count];
        index++;
      }
    }
    data_out = acc;
  }

//------------------------------------------------------------------------------------------------------------------------
// Member Function: firProgCoeffsShiftRegAntiSymmetricEvenTaps()
// Description: firProgCoeffsShiftRegAntiSymmetricEvenTaps() implements a anti-symmetric filter with even number of Taps
// and shift register based implementation

  void firProgCoeffsShiftRegAntiSymmetricEvenTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0.0;
    int ram_addr = 0;
    int index = 0;
    BLK: for ( int i = 0; i < N_TAPS/2; i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH )  {
      index = 0;
#pragma hls_unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
#pragma hls_waive ABR
        acc += (reg[i + index] - reg[N_TAPS - 1 - i - index] ) * coeffs[ram_addr + block_count];
        index++;
      }
    }
    data_out = acc;
  }

//--------------------------------------------------------------------------------------------------------------
// Member Function: firProgCoeffsShiftRegSymmetricOddTaps()
// Description: firProgCoeffsShiftRegSymmetricOddTaps() implements a symmetric filter with odd number of Taps
// and shift register based implementation

  void firProgCoeffsShiftRegSymmetricOddTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0.0;
    ACC_TYPE fold = 0.0;
    int ram_addr = 0;
    int index = 0;
    BLK: for (int i = 0; i < (((N_TAPS - 1) / 2) + 1); i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH) {
      index = 0;
#pragma hls_unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
        if ((i + index)  == (N_TAPS - 1) / 2)
        { fold = ACC_TYPE(reg[i + index]); }                      // Central passes through
        else
        { fold = ACC_TYPE(reg[i + index] + reg[(N_TAPS - 1) - i - index]); }
        index++;
        acc += coeffs[ram_addr + block_count] * fold;
      }
    }
    data_out = acc;
  }

//-----------------------------------------------------------------------------------------------------------------------
// Member Function: firProgCoeffsShiftRegAntiSymmetricOddTaps()
// Description: firProgCoeffsShiftRegAntiSymmetricOddTaps() implements a anti-symmetric filter with odd number of Taps
// and shift register based implementation

  void firProgCoeffsShiftRegAntiSymmetricOddTaps(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
    ACC_TYPE acc = 0.0;
    ACC_TYPE fold = 0.0;
    int ram_addr = 0;
    int index = 0;
    BLK: for (int i = 0; i < (((N_TAPS - 1) / 2) + 1); i+= BLK_SZ, ram_addr+= MEM_WORD_WIDTH) {
      index = 0;
#pragma hls_unroll yes
      MAC: for ( int block_count = BLK_OFFSET; block_count < (BLK_OFFSET + BLK_SZ); block_count ++) {
        if ((i + index)  == (N_TAPS - 1) / 2)
        { fold = ACC_TYPE(reg[i + index]); }                      // Central passes through
        else
        { fold = ACC_TYPE(reg[i + index] - reg[(N_TAPS - 1) - i - index]); }
        index++;
        acc += coeffs[ram_addr + block_count] * fold;
      }
    }
    data_out = acc;
  }

};

#include <mc_scverify.h>

//===================================================================================================================
// Class: ac_fir_reg_share
// Description: This class contains the top function for the register share filter ('run()'). 
//-------------------------------------------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------------------------------
// Member Function: run()
// Description: run() is top function for C++ module.

  // Based on filter type configured it instantiates the desired filter function
  void CCS_BLOCK(run)(IN_TYPE &data_in, COEFF_TYPE coeffs[N_TAPS], OUT_TYPE &data_out) {
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
#pragma hls_waive CNS
    if (ftype == FOLD_ODD) {
      filter.firProgCoeffsShiftRegSymmetricOddTaps(core_in, coeffs, core_out); // Symmetric filter with odd number of Taps
    }
#pragma hls_waive CNS
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

