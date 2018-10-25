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

//***************************************************************************
// File: ac_fft_dif_r2m2p2_sdf.h
//
// Description:
//   Implements an FFT with Mix Radix 2/2^2 with Decimation in Frequency
//   and Single Delay Feedback architecture.
//
// Nomenclature:
//                            ac_fft_dif_r2m2p2_sdf
//                          /    /    |     \       \                               
//                         /   FFT    |   2 mix 2^2  \                               
//                   C-view           |               Single Delay Feedback
//                         Decimation in Frequency
//
//
//  Organization of the design --
//
//  ac_fft_dif_r2m2p2_sdf is a C++ wrapper around a generic C++ Core class 
//  'ac_fft_dif_r2m2p2_sdf_core'. The core class instantiates stages as 
//  'ac_fft_dif_r2m2p2_sdf_stage'. Stages are implemented to
//  behave in Single Delay Feedback structure cascaded in such manner that, it can
//  work for any FFT-Points of 2^n where value of n = 1,2 .. 12.
//  Scaling factor is used to scale down stage output by 1/2 to avoid overflow
//
//   Order of Input/Output
//     Input -- Natural
//     Output-- Bit reversed
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_fft_dif_r2m2p2_sdf.h>
//
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Initialize object of FFT class with 4 FFT points, Memory threshold = 32,
//      // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2.
//      ac_fft_dif_r2m2p2_sdf < 4, 32, 19, 18, 2 > fft_design1;
//      typedef ac_complex<ac_fixed<18, 2, true, AC_TRN, AC_WRAP> > IO_type;
//      // Initialize channels for input and output
//      ac_channel<IO_type> input;
//      ac_channel<IO_type> output;
//      // Make sure to write the inputs to the input port
//      // Call the top-level function.
//      fft_design1.run(input, output);
//
//      CCS_RETURN(0);
//    }
//
// Notes:
//    Attempting to call the function with a type that is not implemented
//    will result in a compile error.
//    Currently, the block only accepts signed ac_complex<ac_fixed> inputs
//    and outputs which use AC_TRN and AC_WRAP as their rounding and
//    overflow modes.
//
//***************************************************************************

#ifndef _INCLUDED_AC_FFT_DIF_R2M2P2_SDF_H_
#define _INCLUDED_AC_FFT_DIF_R2M2P2_SDF_H_

#include <ac_int.h>
#include <ac_fixed.h>
#include <ac_complex.h>

#include <ac_channel.h>
#ifdef ASSERT_ON
#include <ac_assert.h>
#endif
#ifdef COVER_ON
#include <ac_assert.h>
#endif

//=========================================================================
// Class: ac_fft_dif_r2m2p2_sdf_stage
// Description:

template < int STAGE, int MEMORY_THRESHOLD, class com_p, class complex_round, class complex_rnd_ext, class com_mult_type, class complext, class complext_fix >
class ac_fft_dif_r2m2p2_sdf_stage
{
private: // Data
  ac_int < STAGE + 2, false > itrator;
  ac_int < STAGE + 1, false > addr;                        // used for access memory address

  complex_round memory_shift[1 << (STAGE)];                // Memory implementation (must be mapped to RAM)
  complex_round shift_reg[1 << STAGE];                     // Register implementation (must be mapped to registers)

private:

  //----------------------------------------------------------------------
  // Member Function: butterfly
  // Description: Compute Radix-2 Butterfly
  //
  void butterfly(complex_rnd_ext &x, complex_rnd_ext &y) {
    complex_rnd_ext temp_1, temp_2;
    temp_1.r() = x.r();
    temp_1.i() = x.i();
    temp_2.r() = y.r();
    temp_2.i() = y.i();
    x = temp_1 + temp_2;         // Butterfly Computation
    y = temp_1 - temp_2;
  }

  /***************************************************************************************************************************/
  /*
   * function rescale hard-coded to scale down the data by 1/2
   */

  //----------------------------------------------------------------------
  // Member Function: rescale
  // Description:
  //
  complex_round rescale(complex_rnd_ext in) {
    complex_round tx;
    // The function is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    tx.r() = (in.r() >> 1);
    tx.i() = (in.i() >> 1);
    return tx;
  }

  /***************************************************************************************************************************/
  /*
   * function 'multRescale' computes multiplication of twiddle factor and scales by 1/2
   */

  //----------------------------------------------------------------------
  // Member Function: multRescale
  // Description: computes multiplication of twiddle factor and scales
  //   by 1/2.
  //
  complex_round multRescale(complex_rnd_ext x, complext y) {
    com_mult_type temp_1, temp_2, tx;
    complex_round out;
    temp_1.r() = x.r();
    temp_1.i() = x.i();
    temp_2.r() = y.r();
    temp_2.i() = y.i();
    tx = temp_1 * temp_2;
    // The output is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    out.r() = (tx.r() >> 1);
    out.i() = (tx.i() >> 1);
    return out;
  }

  //----------------------------------------------------------------------
  // Member Function: multJ
  // Description: computes multiplication of Complex Input with iota
  //   and Scaling.
  complex_round multJ(complex_rnd_ext x, bool swap) {
    complex_rnd_ext tx;
    complex_round out;

    if (swap) {
      tx.r() = x.i();
      tx.i() = (-x.r());
    } else {
      tx.r() = x.r();
      tx.i() = x.i();
    }

    // The output is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    out.r() = (tx.r() >> 1);
    out.i() = (tx.i() >> 1);

    return out;
  }

  //----------------------------------------------------------------------
  // Member Function: readMem
  // Description: reads data from shift Register or Circular buffer
  //   and Scaling.
  void readMem(complex_rnd_ext &out) {
    if ((1 << STAGE) < (MEMORY_THRESHOLD)) {
      out = shift_reg[(1 << STAGE) - 1];
    } else {
      out = memory_shift[addr & ((1 << STAGE) - 1)];
    }
  }

  //----------------------------------------------------------------------
  // Member Function: writeMem
  // Description: writes data to shift Register or Circular buffer
  //   and Scaling.
  //
  void writeMem(complex_rnd_ext &in) {
    if ((1 << STAGE) < (MEMORY_THRESHOLD)) {
#pragma unroll yes
      for (int i = (1 << STAGE) - 1; i > 0; i--) {
        shift_reg[i] = shift_reg[i - 1];
      }
      shift_reg[0] = in;
    } else {
      memory_shift[addr & ((1 << STAGE) - 1)] = in;
    }
  }

public:

  //----------------------------------------------------------------------
  // Member Function: stageRun
  // Description: compute of stage computation
  //
  void stageRun(com_p &x, com_p &y, const complext w[1 << (STAGE + 1)], bool radix4) {
    complex_rnd_ext xt;
    complex_rnd_ext yt;
    complex_round tmp_out;
    complex_rnd_ext tmp_write_d;
    complext tw, twd;

    ac_int < STAGE + 2, false > n, n_th;
    complext J = complext(0, 1);

    ac_int < STAGE + 1, false > t;
    bool swap, mult;

    xt.r() = x.r();
    xt.i() = x.i();
    readMem(yt);
    mult = 0;

    if (radix4) {               // Switch between Radix 2 and Radix 2^2 butterfly Structure
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
    } else {
      if (itrator[STAGE]) {
        butterfly(yt, xt);
        tmp_out = rescale(yt);
        tmp_write_d = xt;
      } else {
        ac_int < STAGE + 2, false > n;
        n = (itrator) & ((1 << STAGE) - 1);
        complext J = complext(0, 1);

        ac_int < STAGE + 1, false > t;
        // Extraction of twiddle value from 1/8 Cycle of Complex exponential
        t = (1 & ((n << 2) >> STAGE)) ? (ac_int < STAGE + 1, false >) ((((1 << STAGE) >> 2)) - (n & ((((1 << STAGE) >> 2)) - 1))) : (ac_int < STAGE + 1, false >) (n & ((((1 << STAGE) >> 2)) - 1));
        twd = w[t];
        tw.r() = ((1 & ((n << 2) >> STAGE)) | (1 & ((n << 1) >> STAGE))) ? (complext_fix) (-twd.r()) : (complext_fix) (twd.r());
        tw.i() = ((!(1 & ((n << 2) >> STAGE))) | (1 & ((n << 1) >> STAGE))) ? (complext_fix) (twd.i()) : (complext_fix) (-twd.i());
        tw = ((1 & ((n << 2) >> STAGE)) ^ (1 & ((n << 1) >> STAGE))) && (STAGE >= 2) ? (complext) (J * tw.conj()) : (complext) tw;
        tmp_out = multRescale(yt, tw);
        tmp_write_d = xt;
      }
    }
    writeMem(tmp_write_d);
    y.r() = tmp_out.r();
    y.i() = tmp_out.i();
    itrator = itrator + 1;
    addr++;
    return;
  }

  //----------------------------------------------------------------------
  // Constructor
  //   Reset Action will be done here
  //
  ac_fft_dif_r2m2p2_sdf_stage() {
    itrator = 0;
    addr = 0;
    ac::init_array < AC_VAL_DC > (&memory_shift[0], 1 << STAGE);
  }
};

//=========================================================================
// Class: ac_fft_dif_r2m2p2_sdf_core
// Description:
//

template <unsigned N_FFT, int MEMORY_THRESHOLD, int TWID_PREC, int DIF_D0_P, int DIF_D0_I>
class ac_fft_dif_r2m2p2_sdf_core
{
private:
  // Type definitions for Multipliers,Accumulator and stage variable
  typedef ac_fixed<TWID_PREC, 2, true, AC_RND_INF>                     fix_point_tw;
  typedef ac_complex<fix_point_tw>                                     comp_tw;

  typedef ac_fixed<DIF_D0_P, DIF_D0_I, true>                           dif_fix_point;
  typedef ac_complex<dif_fix_point>                                    dif_complex;
  typedef ac_fixed<DIF_D0_P + 1, DIF_D0_I + 1, true, AC_RND, AC_SAT>   dif_fix_round;
  typedef ac_complex<dif_fix_round>                                  	 dif_complex_round;
  typedef ac_fixed<DIF_D0_P + 1, DIF_D0_I + 1, true>                 	 dif_fix_round_ext;
  typedef ac_complex<dif_fix_round_ext>                                dif_complex_round_ext;
  typedef ac_fixed<DIF_D0_P + TWID_PREC - 1, 1 + DIF_D0_I, true>       dif_fix_mul;
  typedef ac_complex<dif_fix_mul>                                      dif_comp_mul;

  // Create stage objects for FFT
  ac_fft_dif_r2m2p2_sdf_stage <11, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_11;
  ac_fft_dif_r2m2p2_sdf_stage <10, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_10;
  ac_fft_dif_r2m2p2_sdf_stage < 9, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_9;
  ac_fft_dif_r2m2p2_sdf_stage < 8, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_8;
  ac_fft_dif_r2m2p2_sdf_stage < 7, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_7;
  ac_fft_dif_r2m2p2_sdf_stage < 6, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_6;
  ac_fft_dif_r2m2p2_sdf_stage < 5, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_5;
  ac_fft_dif_r2m2p2_sdf_stage < 4, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_4;
  ac_fft_dif_r2m2p2_sdf_stage < 3, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_3;
  ac_fft_dif_r2m2p2_sdf_stage < 2, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_2;
  ac_fft_dif_r2m2p2_sdf_stage < 1, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_1;
  ac_fft_dif_r2m2p2_sdf_stage < 0, MEMORY_THRESHOLD, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, comp_tw, fix_point_tw> stage_0;

public:

  //----------------------------------------------------------------------
  // Member Function: fftDifR2m2p2SdfCore
  // Description:
  //   Route data based on number of FFT points and compute stages
  //

  void fftDifR2m2p2SdfCore(const dif_complex &input, dif_complex &output) {
    bool stage_init = false;
    bool r4 = false;

#include <twiddles_64bits.h>

    dif_complex temp_11, temp_10, temp_9, temp_8, temp_7, temp_6, temp_5,
                temp_4, temp_3, temp_2, temp_1, temp_0;

    dif_complex_round temp_rnd_11, temp_rnd_10, temp_rnd_9, temp_rnd_8,
                      temp_rnd_7, temp_rnd_6, temp_rnd_5, temp_rnd_4, temp_rnd_3,
                      temp_rnd_2, temp_rnd_1, temp_rnd_0, out_rnd;

    if (N_FFT >= 4096) {
      temp_rnd_11.r() = input.r();
      temp_rnd_11.i() = input.i();
      temp_11.r() = temp_rnd_11.r();
      temp_11.i() = temp_rnd_11.i();
      stage_init = true;
      r4 = true;
      stage_11.stageRun(temp_11, temp_11, twiddle_11, r4);
    }
    if (N_FFT >= 2048) {
      if (stage_init) {
        temp_rnd_10.r() = temp_11.r();
        temp_rnd_10.i() = temp_11.i();
      } else {
        temp_rnd_10.r() = input.r();
        temp_rnd_10.i() = input.i();
      }
      temp_10.r() = temp_rnd_10.r();
      temp_10.i() = temp_rnd_10.i();
      stage_init = true;
      stage_10.stageRun(temp_10, temp_10, r4 ? twiddle_11 : twiddle_10, r4);
      r4 = true;
    }

    if (N_FFT >= 1024) {
      if (stage_init) {
        temp_rnd_9.r() = temp_10.r();
        temp_rnd_9.i() = temp_10.i();
      } else {
        temp_rnd_9.r() = input.r();
        temp_rnd_9.i() = input.i();
      }
      temp_9.r() = temp_rnd_9.r();
      temp_9.i() = temp_rnd_9.i();
      stage_init = true;
      r4 = true;
      stage_9.stageRun(temp_9, temp_9, twiddle_9, r4);
    }

    if (N_FFT >= 512) {
      if (stage_init) {
        temp_rnd_8.r() = temp_9.r();
        temp_rnd_8.i() = temp_9.i();
      } else {
        temp_rnd_8.r() = input.r();
        temp_rnd_8.i() = input.i();
      }
      temp_8.r() = temp_rnd_8.r();
      temp_8.i() = temp_rnd_8.i();
      stage_init = true;
      stage_8.stageRun(temp_8, temp_8, r4 ? twiddle_9 : twiddle_8, r4);
      r4 = true;
    }

    if (N_FFT >= 256) {
      if (stage_init) {
        temp_rnd_7.r() = temp_8.r();
        temp_rnd_7.i() = temp_8.i();
      } else {
        temp_rnd_7.r() = input.r();
        temp_rnd_7.i() = input.i();
      }
      temp_7.r() = temp_rnd_7.r();
      temp_7.i() = temp_rnd_7.i();
      stage_init = true;
      r4 = true;
      stage_7.stageRun(temp_7, temp_7, twiddle_7,r4);
    }

    if (N_FFT >= 128) {
      if (stage_init) {
        temp_rnd_6.r() = temp_7.r();
        temp_rnd_6.i() = temp_7.i();
      } else {
        temp_rnd_6.r() = input.r();
        temp_rnd_6.i() = input.i();
      }
      temp_6.r() = temp_rnd_6.r();
      temp_6.i() = temp_rnd_6.i();
      stage_init = true;
      stage_6.stageRun(temp_6, temp_6, r4 ? twiddle_7 : twiddle_6, r4);
      r4 = true;
    }

    if (N_FFT >= 64) {
      if (stage_init) {
        temp_rnd_5.r() = temp_6.r();
        temp_rnd_5.i() = temp_6.i();
      } else {
        temp_rnd_5.r() = input.r();
        temp_rnd_5.i() = input.i();
      }
      temp_5.r() = temp_rnd_5.r();
      temp_5.i() = temp_rnd_5.i();
      stage_init = true;
      r4 = true;
      stage_5.stageRun(temp_5, temp_5, twiddle_5, r4);
    }

    if (N_FFT >= 32) {
      if (stage_init) {
        temp_rnd_4.r() = temp_5.r();
        temp_rnd_4.i() = temp_5.i();
      } else {
        temp_rnd_4.r() = input.r();
        temp_rnd_4.i() = input.i();
      }
      temp_4.r() = temp_rnd_4.r();
      temp_4.i() = temp_rnd_4.i();
      stage_init = true;
      stage_4.stageRun(temp_4, temp_4, r4 ? twiddle_5 : twiddle_4, r4);
      r4 = true;
    }

    if (N_FFT >= 16) {
      if (stage_init) {
        temp_rnd_3.r() = temp_4.r();
        temp_rnd_3.i() = temp_4.i();
      } else {
        temp_rnd_3.r() = input.r();
        temp_rnd_3.i() = input.i();
      }
      temp_3.r() = temp_rnd_3.r();
      temp_3.i() = temp_rnd_3.i();
      stage_init = true;
      r4 = true;
      stage_3.stageRun(temp_3, temp_3, twiddle_3, r4);
    }

    if (N_FFT >= 8) {
      if (stage_init) {
        temp_rnd_2.r() = temp_3.r();
        temp_rnd_2.i() = temp_3.i();
      } else {
        temp_rnd_2.r() = input.r();
        temp_rnd_2.i() = input.i();
      }
      temp_2.r() = temp_rnd_2.r();
      temp_2.i() = temp_rnd_2.i();
      stage_init = true;
      stage_2.stageRun(temp_2, temp_2, r4 ? twiddle_3 : twiddle_2, r4);
      r4 = true;
    }

    if (N_FFT >= 4) {
      if (stage_init) {
        temp_rnd_1.r() = temp_2.r();
        temp_rnd_1.i() = temp_2.i();
      } else {
        temp_rnd_1.r() = input.r();
        temp_rnd_1.i() = input.i();
      }
      temp_1.r() = temp_rnd_1.r();
      temp_1.i() = temp_rnd_1.i();
      stage_init = true;
      r4 = true;
      stage_1.stageRun(temp_1, temp_1, twiddle_3, r4);
    }

    if (N_FFT >= 2) {
      if (stage_init) {
        temp_rnd_0.r() = temp_1.r();
        temp_rnd_0.i() = temp_1.i();
      } else {
        temp_rnd_0.r() = input.r();
        temp_rnd_0.i() = input.i();
      }
      temp_0.r() = temp_rnd_0.r();
      temp_0.i() = temp_rnd_0.i();
      stage_init = true;
      stage_0.stageRun(temp_0, temp_0, r4 ? twiddle_1 : twiddle_0, r4);
      r4 = true;
    }

    out_rnd.r() = temp_0.r();
    out_rnd.i() = temp_0.i();
    output.r() = out_rnd.r();
    output.i() = out_rnd.i();
    stage_init = false;
  }
};

//=========================================================================
// Class: ac_fft_dif_r2m2p2_sdf
// Description:
//

template<unsigned N_FFT, int MEMORY_THRESHOLD, int TWID_PREC, int DIF_D0_P, int DIF_D0_I>
class ac_fft_dif_r2m2p2_sdf
{
public:
  bool write_out;

  typedef ac_fixed <DIF_D0_P, DIF_D0_I, true> dif_fx0;
  typedef ac_complex <dif_fx0> dif_cplex0;
  int cnt;

  dif_cplex0 b;

  // Instantiate object of ac_fft_dif_r2m2p2_sdf_core
  ac_fft_dif_r2m2p2_sdf_core <N_FFT, MEMORY_THRESHOLD, TWID_PREC, DIF_D0_P, DIF_D0_I> fft;

  // Constructor
  ac_fft_dif_r2m2p2_sdf() {
    write_out = false;
    b = 0;
  }

  //----------------------------------------------------------------------
  // Member Function: run
  // Description: top-level interface to FFT

#pragma hls_design interface
  void run(ac_channel < dif_cplex0 > &x1, ac_channel < dif_cplex0 > &y1) {
    coverAssert();
    dif_cplex0 a,y;

#pragma hls_pipeline_init_interval 1
    sample_loop: for (int i = 0; i < N_FFT; i++) {
      a = x1.read();
      fft.fftDifR2m2p2SdfCore(a, y);                  /* Calling Core of FFT Design */
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
  // Member Function: coverAssert
  // Description:
  //   Used for basic template assert condition. This helps to validate if
  //   object of this class created in user code has right set of
  //   parameters defined for it. Code will assert during compile time if
  //   incorrect template values are used.
  //
  void coverAssert() {
#ifdef ASSERT_ON
    AC_ASSERT(((N_FFT == 2) || (N_FFT == 4) || (N_FFT == 8) || (N_FFT == 16) || (N_FFT == 32) || (N_FFT == 64) || (N_FFT == 128) ||
               (N_FFT == 256) || (N_FFT == 512) || (N_FFT == 1024) || (N_FFT == 2048) || (N_FFT == 4096)), "N_FFT is not a power of two");
    AC_ASSERT(TWIDDLE_PREC <= 32, "Twiddle bitwidth greater than 32");
#endif
#ifdef COVER_ON
    cover(TWID_PREC <= 5);
#endif
  }
};

#endif

