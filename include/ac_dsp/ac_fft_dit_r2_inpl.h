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
// File: ac_fft_dit_r2_inpl.h
//
// Description:
//    Nomenclature:
//                       ac_fft_dit_r2_inpl
//                     /    /    |   \     \                                         
//                    /   FFT    |  Radix-2 \                                        
//              C-view           |            In-place
//                    Decimation in Time
//
//    Organization:
//     The ac_fft_dit_r2_inpl class serves as a C++ interface to the core class. The
//     ac_fft_dit_r2_inpl_core class is generic and intended to work well even with SystemC
//     implementations in addition to the supplied C++ implementation. The core class
//     instantiates the butterfly as an object of the 'ac_fft_dit_r2_inpl_butterfly' class,
//     handles computation of the FFT Flow Graph and also fetches and writes data.
//     The butterfly class implements a Radix-2 butterfly, and the output of it is
//     hard-coded to scale-down by half.
//
//    Order of Input/Output:
//     Input -- Natural
//     Output-- Bit reversed
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_fft_dit_r2_inpl.h>
//
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Initialize object of FFT class with 4 FFT points,
//      // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2.
//      ac_fft_dit_r2_inpl < 4, 19, 18, 2 > fft_design1;
//      typedef ac_complex<ac_fixed<18, 2, true, AC_TRN, AC_WRAP> > IO_type;
//      // Initialize channels for input and output
//      ac_channel<IO_type> input;
//      ac_channel<IO_type> output;
//      // Write 4 inputs to the input port
//      input.write(IO_type( .125, .25));
//      input.write(IO_type( .375, .50));
//      input.write(IO_type(1.525, .75));
//      input.write(IO_type( .125, .75));
//      // Call the top-level function.
//      fft_design1.run(input, output);
//
//      CCS_RETURN(0);
//    }
//
// Notes:
//    Attempting to call the function with a type that is not implemented will result
//    in a compile error.
//    Currently, the block only accepts signed ac_complex<ac_fixed> inputs and outputs which use AC_TRN
//    and AC_WRAP as their rounding and overflow modes.
//
//*********************************************************************************************************

#ifndef _INCLUDED_AC_FFT_DIT_R2_INPL_H_
#define _INCLUDED_AC_FFT_DIT_R2_INPL_H_

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

#ifndef __SYNTHESIS__
#include <iostream>
using namespace std;
#endif

//*************************************************
//
// class "ac_fft_dit_r2_inpl_butterfly"
//
// class has member functions--
// butterflyR2( ) -- Radix-2 Butterfly Structure
// rescale( )     -- Scale by 1/2
// mult( )        -- Multiply with twiddle factor
// compute( )     -- Compute FFT butterfly call
//
//*************************************************

template < class fix_p, class com_p, class complex_round, class com_rnd_ext, class com_mult, class com_tw >
class ac_fft_dit_r2_inpl_butterfly
{
  typedef ac_fixed < 20, 2, true, AC_RND > fxpt;
  typedef ac_complex < fxpt > comx_tw;

private:

  // function 'butterflyR2' performs Radix-2 Butterfly computation

  void butterflyR2(com_rnd_ext &x, com_rnd_ext &y) {
    com_rnd_ext temp_1, temp_2;
    temp_1 = x;
    temp_2 = y;
    x = temp_1 + temp_2;
    y = temp_1 - temp_2;
  }

  // function rescale scale down the data by 1/2

  complex_round rescale(const com_rnd_ext in) {
    complex_round tx;
    // The addition output is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    tx.r() = (in.r() >> 1);
    tx.i() = (in.i() >> 1);
    return tx;
  }

  // function mult performs twiddle factor multiplication

  com_rnd_ext mult(const com_rnd_ext x, const com_tw y) {
    com_mult temp_1, temp_2, tx;
    com_rnd_ext out;
    temp_1 = x;
    temp_2 = y;
    tx = temp_1 * temp_2;
    out = tx;
    return out;
  }

public:

  void compute(com_p &x, com_p &y, const comx_tw w) {
    com_rnd_ext xt, yt;
    com_tw tw;
    xt = x;
    yt = y;
    tw = (com_tw) w;
    xt = mult(xt, tw);
    butterflyR2(yt, xt);
    x = rescale(xt);
    y = rescale(yt);
  }
};

//**********************************************************************************
// class "ac_fft_dit_r2_inpl_core " implements Fast Fourier Transform (FFT) as
// Decimation
// In Time (DIT) Radix-2 (R2) in In-place (inpl) Architecture ac_fft_dit_r2_inpl_core
// instantiate 1 instance of stages (ac_fft_dit_r2_inpl_butterfly)
// class has member functions--
// swap() -- Swaps Data b/w two banks to avoid Address conflict
//
// fftDitR2InplCore() -- Implements core functionality the FFT.
//
//**********************************************************************************

template < unsigned N_FFT, int TWID_PREC, int DIT_D_P, int DIT_D_I >
class ac_fft_dit_r2_inpl_core
{
private:
  // Type definitions for Multipliers,Accumulator and stage variable for fft
  typedef ac_fixed < TWID_PREC, 2, true, AC_RND_INF > tType;
  typedef ac_complex < tType > cx_tType;

  typedef ac_fixed < DIT_D_P, DIT_D_I, true > dType;
  typedef ac_complex < dType > cx_dType;
  typedef ac_fixed < DIT_D_P, DIT_D_I, true, AC_RND, AC_SAT > dround;
  typedef ac_complex < dround > cx_dround;
  typedef ac_fixed < DIT_D_P + TWID_PREC - 1, 1 + DIT_D_I, true > mType;
  typedef ac_complex < mType > cx_mType;
  typedef ac_fixed < DIT_D_P + 1, DIT_D_I + 1, true > b_dround;
  typedef ac_complex < b_dround > cx_b_dround;

  // Used for swapping Data between two bank based on address Mask and Bank ID
  void swap(const int addr_mask, const int id_bank_1, cx_dType &data_1, cx_dType &data_2) {
    cx_dType data_1_tmp;
    bool swap_d;
    data_1_tmp = data_1;
    swap_d = (addr_mask & id_bank_1);
    data_1 = swap_d ? data_2 : data_1;
    data_2 = swap_d ? data_1_tmp : data_2;
  }

  // FFT computation will be done in 'fftDitR2InplCore'

public:
  void fftDitR2InplCore(cx_dType bank_2[N_FFT / 2], cx_dType bank_1[N_FFT / 2]) {
#include <twiddles_20bits.h>

    const int FFT_STAGES = ac::log2_ceil < N_FFT >::val;

    cx_dType data_1, data_2;
    cx_mType t;

    ac_int < FFT_STAGES + 1, false > n1, n2;
    ac_int < FFT_STAGES + 1, false > idx;
    ac_int < FFT_STAGES, false > addr_mask = 1;
    ac_int < FFT_STAGES, 0 > id_bank_1 = 0;
    ac_int < FFT_STAGES - 1, 0 > id_bank_2 = 0;
    n1 = 0;
    n2 = 1;
    idx = N_FFT;
    cx_tType tw, twd;

#pragma hls_pipeline_init_interval 1
    STAGE_LOOP:for (int i = 0; i < FFT_STAGES; i++) {
      id_bank_1 = 0;
      n1 = n2;
      n2 = n2 + n2;
      idx >>= 1;
      SEGMENT_LOOP:for (int j = 0; j < N_FFT / 2; j++) {
        int k = j;

        BUTTERFLY_LOOP:for (int kk = 0; kk < N_FFT / 2; kk++) {
          id_bank_2 = id_bank_1 ^ (-1 << i);
          data_1 = bank_1[id_bank_1];
          data_2 = bank_2[id_bank_2];
          swap(addr_mask >> 1, id_bank_1, data_1, data_2);
          // Buttererfly Class Object
          ac_fft_dit_r2_inpl_butterfly < dType, cx_dType, cx_dround, cx_b_dround, cx_mType, cx_tType > btrfly;

          ac_int < (FFT_STAGES - 1) + 2, false > n;
          n = (j * idx);
          cx_tType J = comx_twiddle(0, 1);
          // Extraction of twiddle value from 1/8 Cycle of Complex exponential
          ac_int < (FFT_STAGES - 1) + 1, false > t;
          t = (1 & ((n << 2) >> (FFT_STAGES - 1))) ? (ac_int < (FFT_STAGES - 1) + 1, false >) ((((1 << (FFT_STAGES - 1)) >> 2)) - (n & ((((1 << (FFT_STAGES - 1)) >> 2)) - 1))) : (ac_int < (FFT_STAGES - 1) + 1, false >) (n & ((((1 << (FFT_STAGES - 1)) >> 2)) - 1));

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

          btrfly.compute(data_2, data_1, tw);

          swap(addr_mask, id_bank_1, data_1, data_2);
          bank_2[id_bank_2] = data_2;
          bank_1[id_bank_1] = data_1;

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
      addr_mask <<= 1;
    }

  }
};

//**********************************************************************************
//
// class "ac_fft_dit_r2_inpl"
//
// class has member functions--
// void run( )         -- Implements FFT C++-View
// void coverAssert( ) -- Contains basic Assert and cover Conditions
//
//**********************************************************************************

template < unsigned N_FFT, int TWIDDLE_PREC, int DIT_D_P, int DIT_D_I >
class ac_fft_dit_r2_inpl
{
private:

  typedef ac_fixed < DIT_D_P, DIT_D_I, true > dit_fxp_data;
  typedef ac_complex < dit_fxp_data > dit_input, dit_output;

  ac_fft_dit_r2_inpl_core < N_FFT, TWIDDLE_PREC, DIT_D_P, DIT_D_I > fft;
  dit_input bank_1[N_FFT / 2], bank_2[N_FFT / 2];

public:

  ac_fft_dit_r2_inpl() { //Constructor
    ac::init_array<AC_VAL_DC>(&bank_1[0], N_FFT/2);
    ac::init_array<AC_VAL_DC>(&bank_2[0], N_FFT/2);
  }

  // coverAssert() used for basic template assert condition. This helps to validate if object
  // of this class created in user code has right set of parameters defined for it. Code will assert
  // during compile time if incorrect template values are used.

  void coverAssert() {
#ifdef ASSERT_ON
    AC_ASSERT(((N_FFT == 2) || (N_FFT == 4) || (N_FFT == 8) || (N_FFT == 16) || (N_FFT == 32) || (N_FFT == 64) || (N_FFT == 128) ||
               (N_FFT == 256) || (N_FFT == 512) || (N_FFT == 1024) || (N_FFT == 2048) || (N_FFT == 4096)), "N_FFT is not a power of two");
    AC_ASSERT(TWIDDLE_PREC <= 32, "Twiddle factor bitwidth greater than 32");
    AC_ASSERT(TWIDDLE_PREC >= 2,  "Twiddle factor bitwidth lesser than 2");
    AC_ASSERT(DIT_D_P >= DIT_D_I,  "Stage integer width lesser than bitwidth");
#endif
#ifdef COVER_ON
    cover(TWIDDLE_PREC <= 5);
#endif
  }

#pragma hls_design interface
  void run(ac_channel < dit_input > &inst, ac_channel < dit_output > &outst) {
    coverAssert();
#pragma hls_pipeline_init_interval 1
    INPUT_LOOP:for (int i = 0; i < N_FFT; i++) {
      if (i & 1) {
        bank_2[(N_FFT / 2) - (i >> 1) - 1] = inst.read();
      } else {
        bank_1[i >> 1] = inst.read();
      }
    }

    fft.fftDitR2InplCore(bank_2, bank_1);                       /* Calling Core of FFT Design */

#pragma hls_pipeline_init_interval 1
    OUTPUT_LOOP:for (int j = 0; j < N_FFT; j++) {

      if (j < N_FFT / 2) {
        outst.write(bank_1[j]);
      } else {
        outst.write(bank_2[j - (N_FFT / 2)]);
      }
    }
  }
};

#endif

