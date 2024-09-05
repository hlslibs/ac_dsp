/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) DSP Library                                        *
 *                                                                        *
 *  Software Version: 3.6                                                 *
 *                                                                        *
 *  Release Date    : Sun Aug 25 18:24:45 PDT 2024                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.6.0                                               *
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
//*********************************************************************************************************
// File: ac_fft_dit_r2_sdf.h
//
// Description:
//  The ac_fft_dit_r2_sdf class serves as a C++ interface to the core class. The
//  ac_fft_dit_r2_sdf_core class is generic and intended to work well even with SystemC
//  implementations in addition to the supplied C++ implementation. The core class
//  instantiates delay stages as an objects of the 'ac_fft_dit_r2_sdf_stage' class.
//  The stages are implemented as a cascaded Single Feedback structure.
//
// Nomenclature:
//       ac_fft_dit_r2_sdf
//               |   |   |
//               |   |   -------- Single Delay Feedback
//               |   ------------ Radix-2
//               ---------------- Decimation in Time
//
// Order of Input/Output:
//  Input  -- Bit reversed
//  Output -- Natural
//
// Usage:
//  A sample testbench and its implementation look like this:
//  #include <ac_fft_dit_r2_sdf.h>
//
//  #include <mc_scverify.h>
//
//  CCS_MAIN(int arg, char **argc)
//  {
//    // Initialize object of FFT class with 4 FFT points, Memory threshold = 32,
//    // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2.
//    ac_fft_dit_r2_sdf < 4, 32, 19, 18, 2 > fft_design1;
//
//    typedef ac_complex<ac_fixed<18, 2, true, AC_TRN, AC_WRAP> > IO_type;
//
//    // Declare channels for input and output
//    ac_channel<IO_type> input;
//    ac_channel<IO_type> output;
//
//    // Write 4 inputs to the input port
//    input.write(IO_type( .125, .25));
//    input.write(IO_type( .375, .50));
//    input.write(IO_type(1.525, .75));
//    input.write(IO_type( .125, .75));
//
//    // Write 4 zero inputs to the input port to flush out the pipeline
//    input.write(IO_type(0, 0));
//    input.write(IO_type(0, 0));
//    input.write(IO_type(0, 0));
//    input.write(IO_type(0, 0));
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
// Revision History:
//    3.4.3  - dgb - Updated compiler checks to work with MS VS 2019
//    3.3.0  - Added CDesignChecker waivers/fixes for ac_dsp IP blocks.
//             Changes made in general:
//               - CNS violations were waived away.
//               - FXD violations were fixed by changing integer literals to floating point literals.
//               - UMR violation was fixed by using init_array to initialize shift register.
//
//*********************************************************************************************************

#ifndef _INCLUDED_AC_FFT_DIT_R2_SDF_H_
#define _INCLUDED_AC_FFT_DIT_R2_SDF_H_

#include <ac_int.h>
#include <ac_fixed.h>
#include <ac_complex.h>
#include <ac_channel.h>

#ifdef ASSERT_ON
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

#ifndef __SYNTHESIS__
#include <iostream>
using namespace std;
#endif

//=========================================================================
// Class: ac_fft_dit_r2_sdf_stage
// Description: Templatized class for SDF stage computation.
//-------------------------------------------------------------------------

template <int STAGE, int MEM_TH, class com_p, class complex_round, class com_rnd_ext, class com_mult_type, class complext, class complext_fix >
class ac_fft_dit_r2_sdf_stage
{
public:
  //-------------------------------------------------------------------------
  // Member Function: stageRun
  // Description: Top-level function of class, computes output for a stage
  // and calls private member functions.
  //
  void stageRun(com_p &x, const complext w[1 << STAGE]) {
    com_rnd_ext xt;
    com_rnd_ext yt;
    complex_round tmp_out, tmp_res;
    com_rnd_ext tmp_write_d;
    complext tw, twd;

    xt = x;
    readMem(yt);
    if (iterator[STAGE]) {
      ac_int < STAGE + 2, false > n;
      n = (iterator) & ((1 << STAGE) - 1);
      complext J = complext(0.0, 1.0);

      ac_int < STAGE + 1, false > t;
      t = (1 & ((n << 2) >> STAGE)) ? (ac_int < STAGE + 1, false >) ((((1 << STAGE) >> 2)) - (n & ((((1 << STAGE) >> 2)) - 1))) : (ac_int < STAGE + 1, false >) (n & ((((1 << STAGE) >> 2)) - 1));
      twd = w[t];
      tw.r() = ((1 & ((n << 2) >> STAGE)) | (1 & ((n << 1) >> STAGE))) ? (complext_fix) (-twd.r()) : (complext_fix) (twd.r());
      tw.i() = ((!(1 & ((n << 2) >> STAGE))) | (1 & ((n << 1) >> STAGE))) ? (complext_fix) (twd.i()) : (complext_fix) (-twd.i());
      tw = ((1 & ((n << 2) >> STAGE)) ^ (1 & ((n << 1) >> STAGE))) && (STAGE >= 2) ? (complext) (J * tw.conj()) : (complext) tw;
      xt = mult(xt, tw);
      butterfly(yt, xt);
      tmp_res = rescale(xt);
      tmp_write_d = tmp_res;
      tmp_out = rescale(yt);
    } else {
      tmp_out = yt;
      tmp_write_d = xt;
    }

    writeMem(tmp_write_d);

    x = tmp_out;

    iterator++;
    addr++;
    return;
  }

  //-------------------------------------------------------------------------
  // Constructor
  //
  ac_fft_dit_r2_sdf_stage() {
    ac::init_array < AC_VAL_DC > (&memory_shift[0], 1 << STAGE);
    ac::init_array < AC_VAL_DC > (&shift_reg[0], 1 << STAGE);
    // Reset Action here:
    iterator = (1 << STAGE) ^ 1;
    addr = 0;
  }

private:
  ac_int < STAGE + 1, false > iterator;
  ac_int < STAGE + 1, false > addr;             // used for access memory address

  complex_round memory_shift[1 << (STAGE)];     // Memory implementation (must be mapped to RAM)
  complex_round shift_reg[1 << STAGE];          // Register implementation (must be mapped to registers)

  //-------------------------------------------------------------------------
  // Member Function: butterfly
  // Description: Computes Radix-2 Butterfly.
  //
  void butterfly(com_rnd_ext &x, com_rnd_ext &y) {
    com_rnd_ext temp_1, temp_2;
    temp_1 = x;
    temp_2 = y;
    x = temp_1 + temp_2;
    y = temp_1 - temp_2;
  }

  //-------------------------------------------------------------------------
  // Member Function: rescale
  // Description: Function that is hard-coded to scale output data by 1/2.
  //
  complex_round rescale(const com_rnd_ext in) {
    complex_round tx;
    // The function is hard-coded to scale real and imaginary inputs to it by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    tx.r() = (in.r() >> 1);
    tx.i() = (in.i() >> 1);
    return tx;
  }

  //-------------------------------------------------------------------------
  // Member Function: mult
  // Description: Performs twiddle factor multiplication.
  //
  com_rnd_ext mult(const com_rnd_ext x, const complext y) {
    com_mult_type temp_1, temp_2, tx;
    com_rnd_ext out;
    temp_1 = x;
    temp_2 = y;
    tx = temp_1 * temp_2;
    out = tx;
    return out;
  }

  //-------------------------------------------------------------------------
  // Member Function: readMem
  // Description: Reads data from shift register or circular buffer.
  //
  void readMem(com_rnd_ext &out) {
    if ((1 << STAGE) < (MEM_TH)) {
      out = shift_reg[(1 << STAGE) - 1];
    } else {
      out = memory_shift[addr & ((1 << STAGE) - 1)];
    }
  }

  //-------------------------------------------------------------------------
  // Member Function: writeMem
  // Description: Writes data to shift register or circular buffer.
  //
  void writeMem(com_rnd_ext &in) {
    if ((1 << STAGE) < (MEM_TH)) {
#pragma hls_unroll yes
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
// Class: ac_fft_dit_r2_sdf_core
// Description: Instantiates FFT stages and contains core computation for
// FFT design. Can Instantiate upto 12 stages.
//-------------------------------------------------------------------------

template <  unsigned N_FFT, int MEM_TH, int TWID_PREC, int DIT_D0_P, int DIT_D0_I >
class ac_fft_dit_r2_sdf_core
{
private:
  // Typedefs for public function args declared first, to avoid compile-time errors.
  typedef ac_fixed < DIT_D0_P, DIT_D0_I, true > dit_fix_point;
  typedef ac_complex < dit_fix_point >          dit_complex;

public:
  //-------------------------------------------------------------------------
  // Member Function: fftDitR2SdfCore
  // Description: Implements core functionality of FFT. Executes all FFT
  // stages per function call.
  //
  void fftDitR2SdfCore(const dit_complex &input, dit_complex &output) {
#include "twiddles_20bits.h"

    dit_complex temp_11, temp_10, temp_9, temp_8, temp_7, temp_6, temp_5, temp_4, temp_3, temp_2, temp_1, temp_0;

    temp_11 =   (dit_complex_round)input;
#pragma hls_waive CNS
    if (N_FFT >=    2) { stage_11.stageRun(temp_11,  twiddle_0); }

    temp_10 = (dit_complex_round)temp_11;
#pragma hls_waive CNS
    if (N_FFT >=    4) { stage_10.stageRun(temp_10,  twiddle_1); }

    temp_9  = (dit_complex_round)temp_10;
#pragma hls_waive CNS
    if (N_FFT >=    8) {  stage_9.stageRun( temp_9,  twiddle_2); }

    temp_8  =  (dit_complex_round)temp_9;
#pragma hls_waive CNS
    if (N_FFT >=   16) {  stage_8.stageRun( temp_8,  twiddle_3); }

    temp_7  =  (dit_complex_round)temp_8;
#pragma hls_waive CNS
    if (N_FFT >=   32) {  stage_7.stageRun( temp_7,  twiddle_4); }

    temp_6  =  (dit_complex_round)temp_7;
#pragma hls_waive CNS
    if (N_FFT >=   64) {  stage_6.stageRun( temp_6,  twiddle_5); }

    temp_5  =  (dit_complex_round)temp_6;
#pragma hls_waive CNS
    if (N_FFT >=  128) {  stage_5.stageRun( temp_5,  twiddle_6); }

    temp_4  =  (dit_complex_round)temp_5;
#pragma hls_waive CNS
    if (N_FFT >=  256) {  stage_4.stageRun( temp_4,  twiddle_7); }

    temp_3  =  (dit_complex_round)temp_4;
#pragma hls_waive CNS
    if (N_FFT >=  512) {  stage_3.stageRun( temp_3,  twiddle_8); }

    temp_2  =  (dit_complex_round)temp_3;
#pragma hls_waive CNS
    if (N_FFT >= 1024) {  stage_2.stageRun( temp_2,  twiddle_9); }

    temp_1  =  (dit_complex_round)temp_2;
#pragma hls_waive CNS
    if (N_FFT >= 2048) {  stage_1.stageRun( temp_1, twiddle_10); }

    temp_0  =  (dit_complex_round)temp_1;
#pragma hls_waive CNS
    if (N_FFT >= 4096) {  stage_0.stageRun( temp_0, twiddle_11); }

    output = (dit_complex_round)temp_0;
  }

private:
  // Type definitions for Multipliers, Accumulator and stage variable for all stages
  // based on template args
  typedef ac_fixed < TWID_PREC, 2, true, AC_RND_INF >               fix_point_tw;
  typedef ac_complex < fix_point_tw >                               complext;
  typedef ac_fixed < DIT_D0_P, DIT_D0_I, true, AC_RND, AC_SAT >     dit_fix_round;
  typedef ac_complex < dit_fix_round >                              dit_complex_round;
  typedef ac_fixed < DIT_D0_P + 1, DIT_D0_I + 1, true >             dit_fix_round_b;
  typedef ac_complex < dit_fix_round_b >                            dit_complex_round_b;
  typedef ac_fixed < DIT_D0_P + TWID_PREC - 1, 1 + DIT_D0_I, true > dit_fix_mul;
  typedef ac_complex < dit_fix_mul >                                dit_comp_mul;

  // creating stage objects for FFT
  ac_fft_dit_r2_sdf_stage <  0, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_11;
  ac_fft_dit_r2_sdf_stage <  1, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_10;
  ac_fft_dit_r2_sdf_stage <  2, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_9;
  ac_fft_dit_r2_sdf_stage <  3, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_8;
  ac_fft_dit_r2_sdf_stage <  4, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_7;
  ac_fft_dit_r2_sdf_stage <  5, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_6;
  ac_fft_dit_r2_sdf_stage <  6, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_5;
  ac_fft_dit_r2_sdf_stage <  7, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_4;
  ac_fft_dit_r2_sdf_stage <  8, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_3;
  ac_fft_dit_r2_sdf_stage <  9, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_2;
  ac_fft_dit_r2_sdf_stage < 10, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_1;
  ac_fft_dit_r2_sdf_stage < 11, MEM_TH, dit_complex, dit_complex_round, dit_complex_round_b, dit_comp_mul, complext, fix_point_tw > stage_0;
};

//==================================================================================
// Class: ac_fft_dit_r2_sdf
// Description:
// Top-level class of design, instantiated by testbench.
// HLS Interface: run()
//----------------------------------------------------------------------------------

template <unsigned N_FFT, int MEM_TH, int TWID_PREC, int DIT_D0_P, int DIT_D0_I>
class ac_fft_dit_r2_sdf
{
public:
  // Typedefs for public function args declared first, to avoid compile-time errors.
  typedef ac_fixed < DIT_D0_P, DIT_D0_I, true > dit_fx0;
  typedef ac_complex < dit_fx0 > comp_dit0;

  //-------------------------------------------------------------------------
  // Member Function: run
  // Description: Top-level, interface function which is meant to be called
  // by the testbench by instantiating an object of the class with the
  // appropriate template parameters.
  //
#pragma hls_design interface
#pragma hls_pipeline_init_interval 1
  void CCS_BLOCK(run)(ac_channel < comp_dit0 > &x1, ac_channel < comp_dit0 > &y1) {
    coverAssert();
    comp_dit0 a, y;

#pragma hls_pipeline_init_interval 1
    sample_loop: for (int i = 0; i < N_FFT; i++) {
      a = x1.read();
      fft.fftDitR2SdfCore(a, y);                  // Calling Core of FFT Design
      if (write_out) {
        y1.write(b);
      } else {
        b.r() = 0.0;
        b.i() = 0.0;
      }
      b = y;
    }
    write_out = true;
  }

  //-------------------------------------------------------------------------
  // Member Function: coverAssert
  // Description: used for basic template assert condition. This helps to 
  // validate if object of this class created by the user has right set of
  // parameters defined for it. Code will assert during compile-time if 
  // incorrect template values are used.
  //
  void coverAssert() {
#ifdef ASSERT_ON
    static_assert(N_FFT == 2   || N_FFT == 4   || N_FFT == 8    || N_FFT == 16   || N_FFT == 32   || N_FFT == 64 || N_FFT == 128 ||
                  N_FFT == 256 || N_FFT == 512 || N_FFT == 1024 || N_FFT == 2048 || N_FFT == 4096, "N_FFT is not a power of two");
    static_assert(TWID_PREC <= 32, "Twiddle bitwidth greater than 32");
#endif
#ifdef COVER_ON
    cover(TWID_PREC <= 5);
#endif
  }

  //-------------------------------------------------------------------------
  // Constructor
  //
  ac_fft_dit_r2_sdf() {
    // Reset action here:
    write_out = false;
    b.r() = 0.0;
    b.i() = 0.0;
  };

private:
  bool write_out;
  comp_dit0 b;
  ac_fft_dit_r2_sdf_core <N_FFT, MEM_TH, TWID_PREC, DIT_D0_P, DIT_D0_I > fft;
};

#endif

