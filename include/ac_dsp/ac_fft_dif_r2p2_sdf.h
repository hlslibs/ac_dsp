/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) DSP Library                                        *
 *                                                                        *
 *  Software Version: 3.4                                                 *
 *                                                                        *
 *  Release Date    : Sat Jan 23 14:58:27 PST 2021                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.4.0                                               *
 *                                                                        *
 *  Copyright , Mentor Graphics Corporation,                     *
 *                                                                        *
 *  All Rights Reserved.                                                  *
 *  
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
//*********************************************************************************************************
// File:       ac_fft_dif_r2p2_sdf.h
// Description:
//  Implements an FFT with Radix-2^2, Decimation in Frequency
//  and Single Delay Feedback architecture.
//
//  ac_fft_dif_r2p2_sdf is a C++ wrapper around a generic C++ Core class
//  ac_fft_dif_r2p2_sdf_core. The core class instantiates stages as
//  ac_fft_dif_r2p2_sdf_stage. Stages are implemented to
//  behave in Single Delay Feedback structure; cascaded in such manner that they can
//  work for any FFT-Points of 4^n where value of n = 1,2 .. 6.
//  Scaling factor is used to scale down stage output by 1/2 to avoid overflow
//
// Nomenclature:
//       ac_fft_dif_r2p2_sdf
//               |   |   |
//               |   |   -------- Single Delay Feedback
//               |   ------------ Radix-2^2
//               ---------------- Decimation in Frequency
//
// Order of Input/Output:
//  Input  -- Natural
//  Output -- Bit reversed
//
// Usage:
//  A sample testbench and its implementation look like this:
//  #include <ac_fft_dif_r2p2_sdf.h>
//
//  #include <mc_scverify.h>
//
//  CCS_MAIN(int arg, char **argc)
//  {
//    // Initialize object of FFT class with 4 FFT points, Memory threshold = 32,
//    // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2.
//    ac_fft_dif_r2p2_sdf<4, 32, 19, 18, 2> fft_design1;
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
//  Attempting to call the function with a type that is not implemented will result
//  in a compile error.
//  Currently, the block only accepts signed ac_complex<ac_fixed> inputs and outputs which use AC_TRN
//  and AC_WRAP as their rounding and overflow modes.
//
//*********************************************************************************************************

#ifndef _INCLUDED_AC_FFT_DIF_R2P2_SDF_H_
#define _INCLUDED_AC_FFT_DIF_R2P2_SDF_H_

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
#if defined(ASSERT_ON) && __cplusplus < 201103L
#error Please use C++11 or a later standard for compilation.
#endif

#include <mc_scverify.h>

#ifndef __SYNTHESIS__
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Class: ac_fft_dif_r2p2_sdf_stage
// Description: templatized class for a stage
//-------------------------------------------------------------------------

template < int STAGE, int MEM_TH, class com_p, class complex_round, class complex_rnd_ext, class com_mult_type, class complext, class complext_fix >
class ac_fft_dif_r2p2_sdf_stage
{
public:
  //----------------------------------------------------------------------
  // Constructor
  //
  ac_fft_dif_r2p2_sdf_stage() {
    ac::init_array < AC_VAL_DC > (&memory_shift[0], 1 << STAGE);
    // Reset Action Here:
    itrator = 0;
    addr = 0;
  }
  
  //----------------------------------------------------------------------
  // Member Function: stageRun
  // Description: compute of stage computation
  //
  void stageRun(com_p &x, const complext w[1 << (STAGE + 1)]) {
    complex_rnd_ext xt;
    complex_rnd_ext yt;
    complex_round tmp_out;
    complex_rnd_ext tmp_write_d;
    complext tw, twd;

    ac_int < STAGE + 2, false > n, n_th;
    complext J = complext(0, 1);

    ac_int < STAGE + 1, false > t;

    bool swap, mult;

    xt = x;

    readMem(yt);
    mult = 0;
    if (itrator[STAGE]) {

      butterfly(yt, xt);
      tmp_write_d = xt;

      if ((!(STAGE & 1)) && (!itrator[STAGE + 1])) {
        n = (itrator) & (((1 << STAGE)) - 1);
        n_th = n;
        mult = 1;
      } else {
        tmp_out = rescale(yt);
      }
    } else {
      if (!(STAGE & 1)) {
        n = (itrator) & ((1 << STAGE) - 1);
        n_th = itrator[STAGE + 1] ? (ac_int < STAGE + 2, false >) (3 * n) : (ac_int < STAGE + 2, false >) (2 * n);
        mult = 1;
      } else {
        swap = STAGE > 0 ? (bool) itrator[STAGE - 1] : 0;
        tmp_out = multJ(yt, swap);
      }
      tmp_write_d = xt;
    }
    if (mult) {
      t = (1 & ((n_th << 1) >> STAGE)) ? (ac_int < STAGE + 1, false >) ((((1 << (STAGE + 1)) >> 2)) - (n_th & ((((1 << (STAGE + 1)) >> 2)) - 1))) : (ac_int < STAGE + 1, false >) (n_th & ((((1 << (STAGE + 1)) >> 2)) - 1));
      twd = w[t];
      tw.r() = ((1 & ((n_th << 1) >> STAGE)) | (1 & ((n_th) >> STAGE))) ? (complext_fix) (-twd.r()) : (complext_fix) (twd.r());
      tw.i() = ((!(1 & ((n_th << 1) >> STAGE))) | (1 & ((n_th) >> STAGE))) ? (complext_fix) (twd.i()) : (complext_fix) (-twd.i());
      tw = ((1 & ((n_th << 1) >> STAGE)) ^ (1 & ((n_th) >> STAGE))) && ((STAGE + 1) >= 2) ? (complext) (J * tw.conj()) : (complext) tw;
      tw = n_th[STAGE + 1] ? (complext) (-tw) : (complext) (tw);
      tmp_out = multRescale(yt, tw);
    }
    writeMem(tmp_write_d);

    x = tmp_out;
    itrator = itrator + 1;
    addr++;

    return;
  }

private:
  ac_int < STAGE + 2, false > itrator;
  ac_int < STAGE + 1, false > addr;               // used for access memory address
  complex_round memory_shift[1 << (STAGE)];       // Memory implementation (must be mapped to RAM)
  complex_round shift_reg[1 << STAGE];            // Register implementation (must be mapped to registers)
  
  //----------------------------------------------------------------------
  // Member Function: butterfly
  // Description: Performs Radix-2 Butterfly computation
  //
  void butterfly(complex_rnd_ext &x, complex_rnd_ext &y) {
    complex_rnd_ext temp_1, temp_2;

    temp_1 = x;
    temp_2 = y;

    x = temp_1 + temp_2; // Butterfly Computation
    y = temp_1 - temp_2;

  }

  //----------------------------------------------------------------------
  // Member Function: rescale
  // Description: rescale hard-coded to scale down the data by 1/2
  //
  complex_round rescale(complex_rnd_ext in) {
    complex_round tx;

    // The function is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    tx.r() = (in.r() >> 1);
    tx.i() = (in.i() >> 1);

    return tx;
  }

  //----------------------------------------------------------------------
  // Member Function: multRescale
  // Description: Computes multiplication of twiddle factor and scales by 1/2
  //
  complex_round multRescale(complex_rnd_ext x, complext y) {
    com_mult_type temp_1, temp_2, tx;
    complex_round out;

    temp_1 = x;
    temp_2 = y;

    tx = temp_1 * temp_2;

    // The output is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    out.r() = (tx.r() >> 1);
    out.i() = (tx.i() >> 1);

    return out;
  }

  //----------------------------------------------------------------------
  // Member Function: multJ
  // Description: Carries out swapping and scaling.
  //
  complex_round multJ(complex_rnd_ext x, bool swap) {
    complex_rnd_ext tx;
    complex_round out;

    if (swap) {
      tx.r() = x.i();
      tx.i() = (-x.r());
    } else {
      tx = x;
    }

    // The output is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    out.r() = (tx.r() >> 1);
    out.i() = (tx.i() >> 1);

    return out;
  }

  //----------------------------------------------------------------------
  // Member Function: readMem
  // Description: read data from shift register or circular buffer
  //
  void readMem(complex_rnd_ext &out) {
    if ((1 << STAGE) < (MEM_TH)) {
      out = shift_reg[(1 << STAGE) - 1];
    } else {
      out = memory_shift[addr & ((1 << STAGE) - 1)];
    }
  }

  //----------------------------------------------------------------------
  // Member Function: writeMem
  // Description: write data to shift register or circular buffer
  //
  void writeMem(complex_rnd_ext &in) {
    if ((1 << STAGE) < (MEM_TH)) {
#pragma unroll yes
      for (int i = (1 << STAGE) - 1; i > 0; i--) {
        shift_reg[i] = shift_reg[i - 1];
      }
      shift_reg[0] = in;
    } else {
      memory_shift[addr & ((1 << STAGE) - 1)] = in;
    }
  }
};

//=========================================================================
// Class: ac_fft_dif_r2p2_sdf_core
// Description: FFT core class - instantiates stage class to form FFT
//-------------------------------------------------------------------------

template <unsigned N_FFT, int MEM_TH, int TWID_PREC, int DIF_D0_P, int DIF_D0_I>
class ac_fft_dif_r2p2_sdf_core
{

private:
  // Typedefs for public function args declared first, to avoid compile-time errors.
  typedef ac_fixed<DIF_D0_P, DIF_D0_I, true> dif_fix_point;
  typedef ac_complex<dif_fix_point>          dif_complex;

public:
  //----------------------------------------------------------------------
  // Member Function: fftDifR2p2SdfCore
  // Description: Will instantiate all stages of FFT according to the
  // template parameter N_FFT (number of FFT-Points) and pass twiddle
  // vector.
  //
  void fftDifR2p2SdfCore(const dif_complex &input, dif_complex &output) {
    bool stage_init = false;

#include "twiddles_64bits.h"

    dif_complex temp_11, temp_10, temp_9, temp_8, temp_7, temp_6, temp_5, temp_4, temp_3, temp_2, temp_1, temp_0;

    if (N_FFT >= 4096) {
      temp_11 =   (dif_complex_round)input;
      stage_init = true;
      stage_11.stageRun(temp_11, twiddle_11);
    }

    if (N_FFT >= 2048) {
      if (stage_init) { temp_10 = (dif_complex_round)temp_11; }
      else            { temp_10 = (dif_complex_round)input;   }
      stage_init = true;
      stage_10.stageRun(temp_10, twiddle_11);
    }

    if (N_FFT >= 1024) {
      if (stage_init) { temp_9 = (dif_complex_round)temp_10; }
      else            { temp_9 = (dif_complex_round)input;   }
      stage_init = true;
      stage_9.stageRun(temp_9, twiddle_10);
    }

    if (N_FFT >= 512) {
      if (stage_init) { temp_8 = (dif_complex_round)temp_9; }
      else            { temp_8 = (dif_complex_round)input;  }
      stage_init = true;
      stage_8.stageRun(temp_8, twiddle_9);
    }

    if (N_FFT >= 256) {
      if (stage_init) { temp_7 = (dif_complex_round)temp_8; }
      else            { temp_7 = (dif_complex_round)input;  }
      stage_init = true;
      stage_7.stageRun(temp_7, twiddle_8);
    }

    if (N_FFT >= 128) {
      if (stage_init) { temp_6 = (dif_complex_round)temp_7; }
      else            { temp_6 = (dif_complex_round)input;  }
      stage_init = true;
      stage_6.stageRun(temp_6, twiddle_7);
    }

    if (N_FFT >= 64) {
      if (stage_init) { temp_5 = (dif_complex_round)temp_6; }
      else            { temp_5 = (dif_complex_round)input;  }
      stage_init = true;
      stage_5.stageRun(temp_5, twiddle_6);
    }

    if (N_FFT >= 32) {
      if (stage_init) { temp_4 = (dif_complex_round)temp_5; }
      else            { temp_4 = (dif_complex_round)input;  }
      stage_init = true;
      stage_4.stageRun(temp_4, twiddle_5);
    }

    if (N_FFT >= 16) {
      if (stage_init) { temp_3 = (dif_complex_round)temp_4; }
      else            { temp_3 = (dif_complex_round)input;  }
      stage_init = true;
      stage_3.stageRun(temp_3, twiddle_4);
    }

    if (N_FFT >= 8) {
      if (stage_init) { temp_2 = (dif_complex_round)temp_3; }
      else            { temp_2 = (dif_complex_round)input;  }
      stage_init = true;
      stage_2.stageRun(temp_2, twiddle_3);
    }

    if (N_FFT >= 4) {
      if (stage_init) { temp_1 = (dif_complex_round)temp_2; }
      else            { temp_1 = (dif_complex_round)input;  }
      stage_init = true;
      stage_1.stageRun(temp_1, twiddle_2);
    }

    if (N_FFT >= 2) {
      if (stage_init) { temp_0 = (dif_complex_round)temp_1; }
      else            { temp_0 = (dif_complex_round)input;  }
      stage_init = true;
      stage_0.stageRun(temp_0, twiddle_1);
    }

    output = (dif_complex_round)temp_0;

    stage_init = false;
  }
  
private:
  typedef ac_fixed<TWID_PREC, 2, true, AC_RND_INF>                     fix_point_tw;
  typedef ac_complex<fix_point_tw>                                     comp_tw;
  typedef ac_fixed<DIF_D0_P + 1, DIF_D0_I + 1, true, AC_RND, AC_SAT>   dif_fix_round;
  typedef ac_complex<dif_fix_round>                                    dif_complex_round;
  typedef ac_fixed<DIF_D0_P + 1, DIF_D0_I + 1, true, AC_TRN, AC_SAT>   dif_fix_round_ext;
  typedef ac_complex<dif_fix_round_ext>                                dif_complex_round_ext;
  typedef ac_fixed<DIF_D0_P + TWID_PREC - 1, 1 + DIF_D0_I, true>       dif_fix_mul;
  typedef ac_complex<dif_fix_mul>                                      dif_comp_mul;

  // creating stage objects for FFT
  ac_fft_dif_r2p2_sdf_stage <11, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_11;
  ac_fft_dif_r2p2_sdf_stage <10, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_10;
  ac_fft_dif_r2p2_sdf_stage < 9, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_9;
  ac_fft_dif_r2p2_sdf_stage < 8, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_8;
  ac_fft_dif_r2p2_sdf_stage < 7, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_7;
  ac_fft_dif_r2p2_sdf_stage < 6, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_6;
  ac_fft_dif_r2p2_sdf_stage < 5, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_5;
  ac_fft_dif_r2p2_sdf_stage < 4, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_4;
  ac_fft_dif_r2p2_sdf_stage < 3, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_3;
  ac_fft_dif_r2p2_sdf_stage < 2, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_2;
  ac_fft_dif_r2p2_sdf_stage < 1, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_1;
  ac_fft_dif_r2p2_sdf_stage < 0, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_0;
};

//==================================================================================
// Class: ac_fft_dif_r2p2_sdf
// Description: Top-level class, instantiated by testbench.
// HLS Interface: run()
// Note: coverAssert() will be enabled if either macro COVER_ON/ASSERT_ON is defined
//----------------------------------------------------------------------------------

template<unsigned N_FFT, int MEM_TH, int TWID_PREC, int DIF_D0_P, int DIF_D0_I>
class ac_fft_dif_r2p2_sdf
{
public:
  // Typedefs for public function args declared first, to avoid compile-time errors.
  typedef ac_fixed <DIF_D0_P, DIF_D0_I, true> dif_fx0;
  typedef ac_complex <dif_fx0>                dif_cplex0;

  //----------------------------------------------------------------------
  // Member Function: run
  // Description: the top-level interface to the FFT.
  //
#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void CCS_BLOCK(run)(ac_channel < dif_cplex0 > &x1, ac_channel < dif_cplex0 > &y1) {
    coverAssert();
    dif_cplex0 a, y;

#pragma hls_pipeline_init_interval 1
    sample_loop: for (int i = 0; i < N_FFT; i++) {
      a = x1.read();
      fft.fftDifR2p2SdfCore(a, y); // Calling Core of FFT Design
      if (write_out) {
        y1.write(b);
      } else {
        b.r() = 0;
        b.i() = 0;
      }
      b = y;
    }
    write_out = true;
  }

  //----------------------------------------------------------------------
  // Constructor
  //
  ac_fft_dif_r2p2_sdf() {
    // Reset Action Here:
    write_out = false;
    b = 0;
  }
  
  //----------------------------------------------------------------------
  // Member Function: coverAssert
  // Description: used to assert based template conditions to help
  // validate if the instance of this class created in user code has the
  // right set of parameters defined for it. Code will assert during
  // run time if incorrect/inconsistent template parameter values
  // are used.
  //
  void coverAssert() {
#ifdef ASSERT_ON
    static_assert(N_FFT == 4 || N_FFT == 16 || N_FFT == 64 || N_FFT == 256 || N_FFT == 1024 || N_FFT == 4096, "N_FFT is not a power of four");
    static_assert(TWID_PREC <= 32, "Twiddle bitwidth greater than 32");
#endif
#ifdef COVER_ON
    cover(TWID_PREC <= 5);
#endif
  }

private:
  bool             write_out;
  dif_cplex0       b;
  ac_fft_dif_r2p2_sdf_core <N_FFT, MEM_TH, TWID_PREC, DIF_D0_P, DIF_D0_I>  fft;
};

#endif

