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
//***************************************************************************
// File: ac_fft_dif_r2_inpl.h
//
// Description:
//  Implements Decimation-In-Frequency, Radix-2, Inplace Architecture FFT.
//  In-place architecture is implemented by reading and writing data from
//  two memory banks. The butterfly block uses a radix of 2.
//
// Nomenclature:
//       ac_fft_dif_r2_inpl
//               |  |   |
//               |  |   --------- In-place Architecture
//               |  ------------- Radix-2
//               ---------------- Decimation in Frequency
//
// Order of Input/Output:
//  Input  -- Natural
//  Output -- Bit reversed
//
// Usage:
//  A sample testbench and its implementation look like this:
//  #include <ac_fft_dif_r2_inpl.h>
//
//  #include <mc_scverify.h>
//
//  CCS_MAIN(int arg, char **argc)
//  {
//    // Initialize object of FFT class with 4 FFT points,
//    // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2.
//    ac_fft_dif_r2_inpl < 4, 19, 18, 2 > fft_design1;
//
//    typedef ac_complex<ac_fixed<18, 2, true, AC_TRN, AC_WRAP> > IO_type;
//
//    // Declare channels for input and output
//    ac_channel<IO_type> input;
//    ac_channel<IO_type> output;
//
//    // Make sure to write the inputs to the input port
//    for (...) input.write(<value>);
//
//    // Call the top-level function.
//    fft_design1.run(input, output);
//
//    CCS_RETURN(0);
//  }
//
// Notes:
//  Attempting to call the function with a type that is not implemented
//  will result in a compile error.
//  Currently, the block only accepts signed ac_complex<ac_fixed> inputs
//  and outputs which use AC_TRN and AC_WRAP as their rounding and
//  overflow modes.
//
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.4.0  - Added CDesignChecker waivers/fixes for ac_dsp IP blocks.
//             Changes made in general:
//               - CNS violations were waived away.
//               - MXS violations was fixed by changing a variable to an 
//                 unsigned int variable to ensure it is of the same signedness 
//                 to the variable it is assigned to.
//
//
//***************************************************************************

#ifndef _INCLUDED_AC_FFT_DIF_R2_INPL_H_
#define _INCLUDED_AC_FFT_DIF_R2_INPL_H_

#include <ac_int.h>
#include <ac_fixed.h>
#include <ac_complex.h>
#include <ac_channel.h>

#ifdef COVER_ON
#include <ac_assert.h>
#endif

// The coverAssert function uses the static_assert() function, which is only supported by C++11 or later
// compiler standards. Hence, the user should be informed if they have defined COVER_ON but not used
// C++11 or a later compiler standard.
#if defined(ASSERT_ON) && (defined(__GNUC__) && (__cplusplus < 201103L))
#error Please use C++11 or a later standard for compilation.
#endif
#if defined(ASSERT_ON) && (defined(_MSC_VER) && (_MSC_VER < 1920) && !defined(__EDG__))
#error Please use Microsoft VS 2019 or a later standard for compilation.
#endif

#include <mc_scverify.h>

//=========================================================================
// Class: ac_fft_dif_r2_inpl_butterfly
// Description: templatized class for a butterfly.
//-------------------------------------------------------------------------

template < class fix_p, class com_p, class complex_round, class com_rnd_ext, class com_mult, class com_tw >
class ac_fft_dif_r2_inpl_butterfly
{
public:
  //----------------------------------------------------------------------
  // Member Function: compute
  // Description: Compute FFT butterfly call
  //
  void compute(com_p &x, com_p &y, const com_tw w) {
    com_rnd_ext xt, yt;
    complex_round tmp_out;
    com_tw tw;
    xt = x;
    yt = y;
    butterflyR2(yt, xt);
    tw = (com_tw) w;
    tmp_out = mult(xt, tw);
    x = tmp_out;
    y = rescale(yt);
  }

private:
  //----------------------------------------------------------------------
  // Member Function: butterflyR2
  // Description: performs Radix-2 Butterfly computation
  //
  void butterflyR2(com_rnd_ext &x, com_rnd_ext &y) {
    com_rnd_ext temp_1, temp_2;
    temp_1 = x;
    temp_2 = y;
    x = temp_1 + temp_2;
    y = temp_1 - temp_2;
  }

  //----------------------------------------------------------------------
  // Member Function: rescale
  // Description: scale down the data by 1/2
  //
  complex_round rescale(const com_rnd_ext in) {
    complex_round tx;
    // The addition output is hard-coded to scale by 1/2. However, the user can avoid scaling by
    // eliminating the right-shift in the following two equations.
    tx.r() = (in.r() >> 1);
    tx.i() = (in.i() >> 1);
    return tx;
  }

  //----------------------------------------------------------------------
  // Member Function: mult
  // Description: performs twiddle factor multiplication
  //
  complex_round mult(const com_rnd_ext x, const com_tw y) {
    com_mult temp_1, temp_2, tx;
    complex_round out;
    temp_1 = x;
    temp_2 = y;
    tx = temp_1 * temp_2;
    // The output is hard-coded to scale by 1/2. However, the user can avoid scaling by
    // eliminating the right-shift in the following two equations.
    out.r() = (tx.r() >> 1);
    out.i() = (tx.i() >> 1);
    return out;
  }
};

//=========================================================================
// Class: ac_fft_dif_r2_inpl_core
// Description: FFT core class - instantiates butterfly class to form FFT.
//-------------------------------------------------------------------------

template < unsigned N_FFT, int TWID_PREC, int DIF_D_P, int DIF_D_I >
class ac_fft_dif_r2_inpl_core
{
private:
  // Typedefs for public function args declared first, to avoid compile-time errors.
  typedef ac_fixed < DIF_D_P, DIF_D_I, true > dType;
  typedef ac_complex < dType > cx_dType;

public:
  //----------------------------------------------------------------------
  // Member Function: fftDifR2InplCore
  // Description: Implements core functionality the FFT
  //
  void fftDifR2InplCore(cx_dType bank_2[N_FFT / 2], cx_dType bank_1[N_FFT / 2]) {

#include "twiddles_20bits.h"

    const int FFT_STAGES = ac::log2_ceil < N_FFT >::val;

    cx_dType data_l, data_u;
    cx_dround x_tmp, data_tmp;
    cx_mType t;

    ac_int < FFT_STAGES + 1, false > n1, n2;
    ac_int < FFT_STAGES + 1, false > idx;
    ac_int < FFT_STAGES, false > addr_mask = N_FFT >> 1;
    ac_int < FFT_STAGES, 0 > id_bank_1 = 0;
    ac_int < FFT_STAGES - 1, 0 > id_bank_2 = 0;
    n1 = N_FFT;
    n2 = 0;
    idx = 1;
    cx_tType twd, tw;

#pragma hls_pipeline_init_interval 1
    STAGE_LOOP: for (int i = FFT_STAGES - 1; i >= 0; i--) {
      id_bank_1 = 0;
      n2 = n1;
      n1 >>= 1;
      SEGMENT_LOOP: for (int j = 0; j < N_FFT / 2; j++) {
        unsigned int k = j;

        BUTTERFLY_LOOP: for (int kk = 0; kk < N_FFT / 2; kk++) {
          id_bank_2 = id_bank_1 ^ (-1 << i);
          data_l = bank_1[id_bank_1];
          data_u = bank_2[id_bank_2];
          swap(addr_mask, id_bank_1, data_l, data_u);
          // Buttererfly Class Object
          ac_fft_dif_r2_inpl_butterfly < dType, cx_dType, cx_dround, cx_b_dround, cx_mType, cx_tType > btrfly;

          ac_int < (FFT_STAGES - 1) + 2, false > n;
          n = (j * idx);
          cx_tType J = comx_twiddle(0, 1);
          // Extraction of twiddle value from 1/8 Cycle of Complex exponential
          ac_int < (FFT_STAGES - 1) + 1, false > t;
          t = (1 & ((n << 2) >> (FFT_STAGES - 1))) ? (ac_int < (FFT_STAGES - 1) + 1, false >) ((((1 << (FFT_STAGES - 1)) >> 2)) - (n & ((((1 << (FFT_STAGES - 1)) >> 2)) - 1))) : (ac_int < (FFT_STAGES - 1) + 1, false >) (n & ((((1 << (FFT_STAGES - 1)) >> 2)) - 1));

#pragma hls_waive CNS
          // Twiddle ROM selection
          switch (FFT_STAGES) {
            case 1:
              twd = twiddle_0[t];
              break;
            case 2:
              twd = twiddle_1[t];
              break;
            case 3:
              twd = twiddle_2[t];
              break;
            case 4:
              twd = twiddle_3[t];
              break;
            case 5:
              twd = twiddle_4[t];
              break;
            case 6:
              twd = twiddle_5[t];
              break;
            case 7:
              twd = twiddle_6[t];
              break;
            case 8:
              twd = twiddle_7[t];
              break;
            case 9:
              twd = twiddle_8[t];
              break;
            case 10:
              twd = twiddle_9[t];
              break;
            case 11:
              twd = twiddle_10[t];
              break;
            case 12:
              twd = twiddle_11[t];
              break;
          }
          tw.r() = ((1 & ((n << 2) >> (FFT_STAGES - 1))) | (1 & ((n << 1) >> (FFT_STAGES - 1)))) ? (tType) (-twd.r()) : (tType) (twd.r());
          tw.i() = ((!(1 & ((n << 2) >> (FFT_STAGES - 1)))) | (1 & ((n << 1) >> (FFT_STAGES - 1)))) ? (tType) (twd.i()) : (tType) (-twd.i());
          tw = ((1 & ((n << 2) >> (FFT_STAGES - 1))) ^ (1 & ((n << 1) >> (FFT_STAGES - 1)))) && ((FFT_STAGES - 1) >= 2) ? (cx_tType) (J * tw.conj()) : (cx_tType) tw;

          btrfly.compute(data_u, data_l, tw);

          swap(addr_mask >> 1, id_bank_1, data_l, data_u);
          bank_2[id_bank_2] = data_u;
          bank_1[id_bank_1] = data_l;
          k += n2;
          id_bank_1 += (1 << i);

          if (id_bank_1[FFT_STAGES - 1]) {
            id_bank_1[FFT_STAGES - 1] = 0;
            id_bank_1 += 1;
          }
          if (k >= N_FFT - 1) { break; }
        }
        if (j == n1 - 1) { break; }
      }
      addr_mask >>= 1;
      idx <<= 1;
    }
  }

private:
  // Type definitions for Multipliers, Accumulator and stage variable for fft
  typedef ac_fixed < TWID_PREC, 2, true, AC_RND_INF > tType;
  typedef ac_complex < tType > cx_tType;
  typedef ac_fixed < DIF_D_P, DIF_D_I, true, AC_RND, AC_SAT > dround;
  typedef ac_complex < dround > cx_dround;
  typedef ac_fixed < DIF_D_P + TWID_PREC - 1, 1 + DIF_D_I, true > mType;
  typedef ac_complex < mType > cx_mType;
  typedef ac_fixed < DIF_D_P + 1, DIF_D_I + 1, true > b_dround;
  typedef ac_complex < b_dround > cx_b_dround;

  //----------------------------------------------------------------------
  // Member Function: swap
  // Description: Used for swapping Data between two bank based on address
  //   Mask and Bank ID
  //
  void swap(const int addr_mask, const int id_bank_1, cx_dType &data_l, cx_dType &data_u) {
    cx_dType data_l_tmp;
    bool swap_d;
    data_l_tmp = data_l;
    swap_d = (addr_mask & id_bank_1);
    data_l = swap_d ? data_u : data_l;
    data_u = swap_d ? data_l_tmp : data_u;
  }
};

//==================================================================================
// Class: ac_fft_dif_r2_inpl
// Description: Top-level class, instantiated by testbench.
// HLS Interface: run()
// Note: coverAssert() will be enabled if either macro COVER_ON/ASSERT_ON is defined
//----------------------------------------------------------------------------------

template < unsigned N_FFT, int TWIDDLE_PREC, int DIF_D_P, int DIF_D_I >
class ac_fft_dif_r2_inpl
{
public:
  // Typedefs for public function args declared first, to avoid compile-time errors.
  typedef ac_fixed<DIF_D_P, DIF_D_I, true> dif_fxp_data;
  typedef ac_complex<dif_fxp_data> dif_input, dif_output;

  //----------------------------------------------------------------------
  // Member Function: run
  // Description: top-level interface to FFT
  //
#pragma hls_design interface
  void CCS_BLOCK(run)(ac_channel<dif_input> &inst, ac_channel<dif_output> &outst) {
    void coverAssert();
#pragma hls_pipeline_init_interval 1
    INPUT_LOOP: for (int j = 0; j < N_FFT; j++) {
      if (j < N_FFT / 2) {
        bank_1[j] = inst.read();
      } else {
        bank_2[j - (N_FFT / 2)] = inst.read();
      }
    }

    fft.fftDifR2InplCore(bank_2, bank_1);  // Calling Core of FFT Design

#pragma hls_pipeline_init_interval 1
    OUTPUT_LOOP: for (int i = 0; i < N_FFT; i++) {
      if (i & 1) {
        outst.write(bank_2[(N_FFT / 2) - (i >> 1) - 1]);
      } else {
        outst.write(bank_1[i >> 1]);
      }
    }
  }

  //----------------------------------------------------------------------
  // Constructor
  //
  ac_fft_dif_r2_inpl() {
    ac::init_array<AC_VAL_DC>(&bank_1[0], N_FFT/2);
    ac::init_array<AC_VAL_DC>(&bank_2[0], N_FFT/2);
  }

  //----------------------------------------------------------------------
  // Member Function: coverAssert
  // Description:
  //   Used for basic template assert condition. This helps to validate if
  //   object of this class created in user code has right set of
  //   parameters defined for it. Code will assert during compile-time if
  //   incorrect template values are used.
  //
  void coverAssert() {
#ifdef ASSERT_ON
    static_assert(N_FFT == 2   || N_FFT == 4   || N_FFT == 8    || N_FFT == 16   || N_FFT == 32   || N_FFT == 64 || N_FFT == 128 ||
                  N_FFT == 256 || N_FFT == 512 || N_FFT == 1024 || N_FFT == 2048 || N_FFT == 4096, "N_FFT is not a power of two");
    static_assert(TWIDDLE_PREC <= 32, "Twiddle bitwidth greater than 32");
    static_assert(TWIDDLE_PREC >= 2,  "Twiddle bitwidth lesser than 2");
    static_assert(DIF_D_P >= DIF_D_I,  "Stage integer width lesser than bitwidth");
#endif
#ifdef COVER_ON
    cover(TWIDDLE_PREC <= 19);
#endif
  }

private:
  // Instantiate object of ac_fft_dif_r2_inpl_core class
  ac_fft_dif_r2_inpl_core < N_FFT, TWIDDLE_PREC, DIF_D_P, DIF_D_I > fft;
  // Allocate banks
  dif_input bank_1[N_FFT / 2], bank_2[N_FFT / 2];
};

#endif

