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

/************************************************************************************/

// File:       ac_fft_dif_r2pX_bfp_inpl.h
// Description:                                 Block Floating Point
// Nomenclature:                               /
//                           ac_fft_dif_r2pX_bfp_inpl
//                         /    /    |   \         \                                
//                        /   FFT    |  R-2^x       \                               
//                  C-view           |               In-place
//                        Decimation in Frequency
//
//
//  Organization of the design --
//
// ac_fft_dif_r2pX_bfp_inpl is a C++  wrapper around a generic Core class 
// 'ac_fft_dif_r2pX_bfp_inpl_core'.  'ac_fft_dif_r2pX_bfp_inpl_core' instantiate butterfly as
// 'ac_fft_dif_r2pX_bfp_inpl_butterfly' and Compute FFT Flow Graph fetching and writing
//  Data from Radix No. of Memory Banks
//
// Architecture view below is showing input and output interface with memory and
// Butterfly Engine output will came out along with a scaling factoe of FFT
//
//
// 'ac_fft_dif_r2pX_bfp_inpl_butterfly'
// implements  Radix-engine comprise. Butterfly output can  scale-down by 1/2 using
// bit pattern of stage_scaling parameter scaling factor provided to butterfly is
// is onthe basis of dynaminc Rgange Analysis done in previous stage
//
//
// Order of Input/Output
// Input -- Natural
// Output-- ORDER = 0 Bit reversed
//                = 1 Natural
//
//
/************************************************************************************/

/************************************************************************************/

/*     +--------+          +--------------------------+          +--------+
       |        |          |                          |          |        | Output
Input  | INPUT  |          |                          |          | OUTPUT +--------->
-------> UNIT   | +-------->          RADIX N         +--------+ |  UNIT  | Stream
Stream |        | |        |        BUTTERFLY         |        | |        | Scaling
       |        | |  +----->          ENGINE          +-----+  | |        +--------->
       |        | |  |     |                          |     |  | |        | Factor
       +---+----+ |  |  +-->                          +--+  |  | +---^----+
           |      |  |  |  |                          |  |  |  |     |
           |      |  |  |  +--------------------------+  |  |  |     |
           |      |  |  |                                |  |  |     |
           |      |  |  |                                |  |  |     |
           |      |  |  |                                |  |  |     |
           |      |  |  |   +-------------------------+  |  |  |     |
           |      |  |  +---+         BANK 1          <--+  |  |     |
           |      |  |      |                         |     |  |     |
           |      |  |      +-------------------------+     |  |     |
           |      |  |      +-------------------------+     |  |     |
           |      |  +------+         BANK 2          <-----+  |     |
           |      |         |                         |        |     |
           |      |         +-------------------------+        |     |
           |      |         +-------------------------+        |     |
           |      +---------+         BANK N          <--------+     |
           |                |                         |              |
           |                +------^-----------+------+              |
           |                       |           |                     |
           +-----------------------+           +---------------------++             */

//*******************************************************************************************************
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_fft_dif_r2pX_bfp_inpl.h>
//
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      typedef ac_complex<ac_fixed<18, 2, true, AC_TRN, AC_WRAP> > I_type;
//      typedef ac_int < ac::log2_ceil < ac::log2_ceil < F_N >::val >::val +1, 1 > scale_typ;
//      struct out_struct {
//      dif_out manti;
//      scale_typ expo;
//      };
//      // Initialize object of FFT class with 4 FFT points, Radix = 2, bit-reversed output,
//      // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2 and the struct out_struct
//      ac_fft_dif_r2pX_dyn_inpl < 4, 2, 0, 19, 18, 2, out_struct > fft_design1;
//      // Initialize channels for input and output
//      ac_channel<I_type> input;
//      ac_channel<out_struct> output;
//      // Make sure to write the inputs to the input port
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


#ifndef _INCLUDED_AC_FFT_DIF_R2PX_BFP_INPL_H_
#define _INCLUDED_AC_FFT_DIF_R2PX_BFP_INPL_H_

#include <ac_int.h>
#include <ac_fixed.h>
#include <ac_complex.h>
#include <ac_channel.h>

#include<iostream>
using namespace std;

#ifdef ASSERT_ON
#include <ac_assert.h>
#endif
#ifdef COVER_ON
#include <ac_assert.h>
#endif

template < int N_FFT, int RADX, int DIF_D_P, int DIF_D_I, class fix_p, class com_p, class complex_round, class rnd_ext, class com_rnd_ext, class com_mult, class com_tw > 
class ac_fft_dif_r2pX_bfp_inpl_butterfly
{
private:
  /*
   * function 'fft_dif_r2pX_bfp_inpl_radix_butterfly' performs Radix Engine performs computation as Radix-2 flow graph
   */

  void fft_dif_r2pX_bfp_inpl_radix_butterfly (int stage_n, com_rnd_ext x[RADX], com_rnd_ext y[RADX], ac_int < ac::log2_ceil < RADX >::val, 0 > scale_fac) {
    const int logN = ac::log2_ceil < N_FFT >::val;
    const int logRad = ac::log2_ceil < RADX >::val;
    const int loglogRad = ac::log2_ceil < logRad >::val;
    const int mixR = logN % logRad;
    com_tw Rw;

#include <twiddlesR_64bits.h>

#pragma unroll yes
    RADIX_STAGE_LOOP: for (ac_int < logRad, 0 > Rstage = logRad - 1; Rstage >= 0; Rstage--) {
      // Radix stage
      com_rnd_ext a[RADX];
#pragma unroll yes
      RESCALE_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
        com_rnd_ext temp;
        bool ch_scale = 1 & (scale_fac >> (logRad - 1 - Rstage));
        if (ch_scale) {
          temp.r () = (x[rad_itr].r () >> 1);
          temp.i () = (x[rad_itr].i () >> 1);
        } else {
          temp.r () = (x[rad_itr].r ());
          temp.i () = (x[rad_itr].i ());
        }
        x[rad_itr].r () = temp.r ();
        x[rad_itr].i () = temp.i ();
      }

#pragma unroll yes
      INPUT_REG_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
        a[rad_itr].r () = x[rad_itr].r ();
        a[rad_itr].i () = x[rad_itr].i ();
      }

#pragma unroll yes
      RADIX_BUTTERFLY_LOOP: for (ac_int < logRad, 0 > B_count = 0; B_count < RADX / 2; B_count++) {
        // Radix butterfly
        ac_int < logRad, 0 > addrs1, addrs2, twid_add;
        twid_add = ((1 << (logRad - 1)) - 1) & (B_count << (logRad - Rstage - 1));
        switch (RADX) {
          case 2:
            Rw = Rtw0[twid_add];
            break;
          case 4:
            Rw = Rtw1[twid_add];
            break;
          case 8:
            Rw = Rtw2[twid_add];
            break;
          case 16:
            Rw = Rtw3[twid_add];
            break;
          case 32:
            Rw = Rtw4[twid_add];
            break;
          case 64:
            Rw = Rtw5[twid_add];
            break;
          case 128:
            Rw = Rtw6[twid_add];
            break;
        }

#pragma unroll yes
        RADIX_ADD_GEN_LOOP: for (int bit_count = 0; bit_count < logRad; bit_count++) {
          if (bit_count < Rstage) {
            addrs1[bit_count] = B_count[bit_count];
            addrs2[bit_count] = B_count[bit_count];
          }
          if (bit_count == Rstage) {
            addrs1[bit_count] = 0;
            addrs2[bit_count] = 1;
          }
          if (bit_count > Rstage) {
            addrs1[bit_count] = B_count[bit_count - 1];
            addrs2[bit_count] = B_count[bit_count - 1];
          }
        }

        if ((((stage_n != 0)) || (!(Rstage >= mixR)) || (mixR == 0))) {
          x[addrs1] = a[addrs1] + a[addrs2];
          x[addrs2] = a[addrs1] - a[addrs2];
        } else {
          x[addrs1] = a[addrs1];
          x[addrs2] = a[addrs2];
        }

        com_mult temp_x, temp_tw;

        temp_x = x[addrs2];
        if ((((stage_n != 0)) || (!(Rstage >= mixR)) || (mixR == 0))) {
          temp_tw = Rw;
        } else {
          temp_tw.r () = 1;
          temp_tw.i () = 0;
        }
        x[addrs2] = (com_rnd_ext) (temp_x * temp_tw);
        x[addrs1] = (com_rnd_ext) (x[addrs1]);
      }

      if (Rstage == 0) { break; }
    }
#pragma unroll yes
    OUTPUT_REG_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      y[rad_itr].r () = x[rad_itr].r ();
      y[rad_itr].i () = x[rad_itr].i ();
    }
  }

  void multRescale (com_rnd_ext x[RADX], com_rnd_ext yo[RADX], com_tw y[RADX], ac_int < ac::log2_ceil < RADX >::val, 0 > scale_fac) {
    const int logRad = ac::log2_ceil < RADX >::val;
    com_mult a[RADX], b[RADX], tx[RADX];
    complex_round temp;
    com_rnd_ext temps;
    bool ch_scale = 1 & (scale_fac >> (logRad - 1));

#pragma unroll yes
    TWIDDLW_MULT_LOOP: for (int rad_itr = 1; rad_itr < RADX; rad_itr++) {

      a[rad_itr].r () = x[rad_itr].r ();
      a[rad_itr].i () = x[rad_itr].i ();
      b[rad_itr].r () = y[rad_itr].r ();
      b[rad_itr].i () = y[rad_itr].i ();

      tx[rad_itr] = a[rad_itr] * b[rad_itr];

      temp.r () = (tx[rad_itr].r ());
      temp.i () = (tx[rad_itr].i ());

      yo[rad_itr].r () = temp.r ();
      yo[rad_itr].i () = temp.i ();
    }

    temps = x[0];
    temp.r () = temps.r ();
    temp.i () = temps.i ();
    yo[0].r () = temp.r ();
    yo[0].i () = temp.i ();
  }

  ac_int < ac::log2_ceil < RADX >::val, 0 > scale_count;

public:

  void compute (int stage_n, com_p x[RADX], com_p y[RADX], const com_tw w[RADX],  ac_int < ac::log2_ceil < RADX >::val, 0 >  &bscale_in) {
    const int logRad = ac::log2_ceil < RADX >::val;
    com_rnd_ext xt[RADX];
    com_rnd_ext xyt[RADX];
    com_rnd_ext yt[RADX];
    complex_round yr[RADX];
    complex_round tmp_out;
    com_tw tw[RADX];
    ac_int < logRad, 0 > scale_fac;
    scale_fac=bscale_in;
#pragma unroll yes
    TWIDDLE_ROUNDING_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      xt[rad_itr].r () = x[rad_itr].r ();
      xt[rad_itr].i () = x[rad_itr].i ();
      tw[rad_itr] = (com_tw) w[rad_itr];
    }

    fft_dif_r2pX_bfp_inpl_radix_butterfly (stage_n, xt, xyt, scale_fac);
    multRescale (xyt, yt, tw, scale_fac);
#pragma unroll yes
    ROUNDING_LOOP_: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      yr[rad_itr].r () = yt[rad_itr].r ();
      yr[rad_itr].i () = yt[rad_itr].i ();
    }

#pragma unroll yes
    OUTPUT_ASSINGMENT_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      y[rad_itr].r () = yr[rad_itr].r ();
      y[rad_itr].i () = yr[rad_itr].i ();
    }
  }
};

template < int N_FFT, int RADIX, int ORDER, int TWID_PREC, int DIF_D_P, int DIF_D_I >
class ac_fft_dif_r2pX_bfp_inpl_core
{
  /*
   * Type definitions for Multipliers,Accumulator and stage variable for fft
   */

public:

  typedef ac_fixed < TWID_PREC, 2, true, AC_RND_INF > tType;
  typedef ac_complex < tType > cx_tType;
  typedef ac_fixed < DIF_D_P, DIF_D_I, true > dType;
  typedef ac_complex < dType > cx_dType;
  typedef ac_fixed < DIF_D_P, DIF_D_I, true, AC_RND, AC_SAT > dround;
  typedef ac_complex < dround > cx_dround;
  typedef ac_fixed < DIF_D_P + TWID_PREC - 1, 1 + DIF_D_I, true > mType;
  typedef ac_complex < mType > cx_mType;
  typedef ac_fixed < DIF_D_P + 1, DIF_D_I + 1, true, AC_RND, AC_SAT > b_dround;
  typedef ac_complex < b_dround > cx_b_dround;
  template < int a > ac_int < a, 0 > bitrevint (ac_int < a, 0 > &Num) {
    ac_int < a, 0 > Num_br;

#pragma unroll yes
    BITREVERSAL_CORE: for (ac_int < a + 1, 0 > itrator = 0; itrator < a; itrator++) {
      Num_br[a - 1 - itrator] = Num[itrator]; // bitreversal ac_int
    }
    return Num_br;
  }
  ac_int<ac::log2_ceil < RADIX >::val, 0> dynRangeCompute (cx_dType butffinp[RADIX]) {
    const int logRx = ac::log2_ceil < RADIX >::val;
    ac_int<logRx, 0> abs_dynRange=0, rangDyn[RADIX];
#pragma unroll yes
    for (int itm = 0; itm < RADIX ; itm++ ) {
      rangDyn[itm]=0;
    }
#pragma unroll yes
    for (int itra = 0; itra < RADIX ; itra++ ) {
      rangDyn[itra]= rangDyn[itra]|(butffinp[itra].r()>=0?butffinp[itra].r().template slc < logRx >(DIF_D_P - logRx - 1 ):(~(butffinp[itra].r().template slc < logRx >(DIF_D_P - logRx - 1 ))))|(butffinp[itra].i()>=0?butffinp[itra].i().template slc < logRx >(DIF_D_P - logRx - 1 ):(~(butffinp[itra].i().template slc < logRx >(DIF_D_P - logRx - 1 ))));
    }
#pragma unroll yes
    for (int itrs = 0; itrs < RADIX ; itrs++ ) {
      abs_dynRange = abs_dynRange|rangDyn[itrs];
    }
    return abs_dynRange;
  }

  void fftDifR2pXInplCore (cx_dType bank[RADIX][N_FFT / RADIX],ac_int <  ac::log2_ceil < DIF_D_P >::val, 1 >  &upscl, ac_int < ac::log2_ceil < ac::log2_ceil < N_FFT >::val >::val +1, 1 > &scale_out ) {
    /*
     * FFT computation will be done in 'fftDifR2pXInplCore'
     */
#include <twiddles_64bits.h>

    const int logN_a = ac::log2_ceil < N_FFT >::val;
    const int loglogN_a = ac::log2_ceil < logN_a >::val;
    const int logRad = ac::log2_ceil < RADIX >::val;
    const int mixR = logN_a % logRad;
    const int RADIX_R = (logRad - mixR) % logRad;
    const int Nstage_a = logN_a / logRad;
    const int modi_sh = (Nstage_a + (mixR != 0)) * logRad;
    const int N_FFT_m = 1 << (modi_sh);
    const int logN = ac::log2_ceil < N_FFT_m >::val;
    const int Nstage = logN / logRad;
    const dType VAL_MAX = (~(dType)(1<<(DIF_D_I-1)));
    cx_dType data_in[RADIX],data_out[RADIX];

    ac_int < logN + 1, false > n1, n2;
    ac_int < logN + 1, false > idx;
    ac_int < logN - logRad + 1, 0 > bank_add_gen = 0;
    ac_int < logN, 0 > xadd = 0;
    ac_int < logN - logRad + 1, 0 > id_bank[RADIX], cur_id_bank, fet_id_data;
    ac_int < loglogN_a, 0> scaling_out ;
    ac_int < logRad, 0 > scaling_fac ;
    ac_int < logN - logRad + 1, 0 > prev_dn_addr = 0, next_dn_addr = 0;
    const ac_int < logN - logRad + 2, 0 > fix_zero_wire_wdth = ((1 << (logN - logRad)) - 1);
    ac_int < ac::log2_ceil < ac::log2_ceil < N_FFT >::val >::val, 0 > dnscl=logRad ;

    n1 = N_FFT_m;
    n2 = 0;
    idx = 1;
    cx_tType twd, tw, tw_vac[RADIX];

    dType  dynRange_temp[RADIX];
    ac_int<ac::log2_ceil < RADIX >::val, 0>  dynRange;
    /*
     * Radix -X engine instantiation
     */
    ac_fft_dif_r2pX_bfp_inpl_butterfly < N_FFT, RADIX, DIF_D_P, DIF_D_I, dType, cx_dType, cx_dround, b_dround, cx_b_dround, cx_mType, cx_tType > btrfly;

    scaling_fac = RADIX-1;

#ifndef __SYNTHESIS__
    std::cout<<"Stage Scale Dynamic Range filled"<<endl;
#endif

    STAGE_LOOP: for (int i = Nstage - 1; i >= 0; i--) {
      // stage
      bank_add_gen = 0;
      n2 = n1;
      n1 >>= logRad;
      ac_int<logRad, 0> abs_dynRange=0;
#pragma unroll yes
      dynRange=0;
#pragma hls_pipeline_init_interval 1
      SEGMENT_LOOP: for (int j = 0; j < N_FFT / RADIX; j++) {
        // segment
        int k = j;

        BUTTERFLY_LOOP: for (int kk = 0; kk < N_FFT / RADIX; kk++) {
          // butterfly

#pragma unroll yes
          BANK_ADD_GENERATOR: for (ac_int < logRad + 1, 0 > m = 0; m < RADIX; m++) {
            ac_int < logRad, 0 > m_no_msb = 0;
            xadd = 0;
            m_no_msb = m;
#pragma unroll yes
            BANK_ADD_GEN_NEST: for (int slc_j = logN - logRad; slc_j >= (i * logRad); slc_j = (slc_j - logRad)) {
              xadd.set_slc (slc_j, m_no_msb);
            }
            id_bank[m] = fix_zero_wire_wdth & (bank_add_gen ^ xadd);
          }

#pragma unroll yes
          DATA_FETCH_FROM_BANKS: for (int mn = 0; mn < RADIX; mn++) {
            int bank_add_fet = mn ^ (bank_add_gen.template slc < logRad > (i * logRad));

            fet_id_data = id_bank[mn] & ((N_FFT / RADIX) - 1);

            data_in[bank_add_fet] = bank[mn][fet_id_data];
            if (i == (Nstage - 1) ) {
              if (upscl>=0) {
                data_in[bank_add_fet].r()= data_in[bank_add_fet].r()<<upscl;
                data_in[bank_add_fet].i()= data_in[bank_add_fet].i()<<upscl;
              } else {
                dnscl=-upscl;
                scaling_fac=(1<<dnscl)-1;
              }
            }
          }
          ac_int < (logN + 1), false > n, n_vac[RADIX];
          n = (j * idx);

#pragma unroll yes
          BITREVERCE_TWIDDLE_ADDRESS: for (ac_int < logRad + 1, 0 > n_itr = 0; n_itr < RADIX; n_itr++) {
            ac_int < logRad, 0 > n_nmsb = n_itr;
            n_nmsb = bitrevint < logRad > (n_nmsb);
            if ((i == Nstage - 1) && (RADIX_R != 0)) {
              n_vac[n_nmsb] = (n_itr >> RADIX_R) * (n + ((n_nmsb >> (logRad - RADIX_R)) * (N_FFT / RADIX)));
            } else {
              n_vac[n_nmsb] = (n_itr * n) >> RADIX_R;
            }

          }
          cx_tType J = comx_twiddle (0, 1);

          tw_vac[0] = comx_twiddle (1, 0);

#pragma unroll yes
          TWIDDLE_VAC_GEN: for (ac_int < logRad + 1, 0 > tw_itr = 1; tw_itr < RADIX; tw_itr++) {
            ac_int < logN, false > t;
            n = n_vac[tw_itr];
            t = (1 & ((n << 2) >> (logN_a - 1))) ? (ac_int < (logN_a - 1) + 1, false >) ((((1 << (logN_a - 1)) >> 2)) - (n & ((((1 << (logN_a - 1)) >> 2)) - 1))) : (ac_int < (logN_a - 1) + 1, false >) (n & ((((1 << (logN_a - 1)) >> 2)) - 1));

            switch (logN_a) {
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
              case 13:
                twd = twiddle_12[t];
                break;
            }
            tw.r () = ((1 & ((n << 2) >> (logN_a - 1))) | (1 & ((n << 1) >> (logN_a - 1)))) ? (tType) (-twd.r ()) : (tType) (twd.r ());
            tw.i () = ((!(1 & ((n << 2) >> (logN_a - 1)))) | (1 & ((n << 1) >> (logN_a - 1)))) ? (tType) (twd.i ()) : (tType) (-twd.i ());
            tw = ((1 & ((n << 2) >> (logN_a - 1))) ^ (1 & ((n << 1) >> (logN_a - 1)))) && ((logN_a - 1) >= 2) ? (cx_tType) (J * tw.conj ()) : (cx_tType) tw;
            tw = n[logN_a - 1] ? (cx_tType) (-tw) : (cx_tType) (tw);
            tw_vac[tw_itr] = tw;
          }
          int temp_count = 0;

          btrfly.compute (i - Nstage + 1, data_in, data_out, tw_vac, scaling_fac );
          /*
           * dynRangeCompute(<data out from butterfly>) provides dynamic logRad MSBs of OR with all real and imaginary Part
           * dynRange of all butterfly is OR to get absolute Range of one Stage to generate scaling factor for next stage.
           */
          dynRange = dynRangeCompute(data_out);
          abs_dynRange = abs_dynRange|dynRange;
#pragma unroll yes
          PUTTING_DATA_IN_BANKS: for (int amn = 0; amn < RADIX; amn++) {
            cur_id_bank = id_bank[amn] & ((N_FFT / RADIX) - 1);
            int data_add = amn ^ (bank_add_gen.template slc < logRad > (i > 0 ? (i - 1) * logRad : logN));

            bank[amn][cur_id_bank] = data_out[data_add];
          }
          k += n2;
          bank_add_gen += (1 << ((i) * logRad));

          if ((bank_add_gen >= N_FFT / RADIX)) {
            // if idx overflows, wrap
            bank_add_gen = bank_add_gen & ((N_FFT / RADIX) - 1);
            bank_add_gen += 1;
          }

          if (k >= N_FFT - 1) { break; }
        }
        if (j == n1 - 1) { break; }
      }

      ac_int<logRad,0> dyn_chk = abs_dynRange;
      ac_int<logRad+1,0> check_d;
      check_d=0;
#pragma unroll yes
      for (int chk = 0; chk < logRad ; chk++ ) {
        check_d+=dyn_chk[chk];
      }
      scaling_fac=dyn_chk;
      if ( i!=0 ) {
        dnscl=check_d+dnscl;
      }
#ifndef __SYNTHESIS__
      std::cout<<"  "<<(Nstage-i)<<"     "<<check_d<<"    "<<((100*(DIF_D_P-logRad+check_d)/DIF_D_P))<<"%"<<endl;
#endif
      idx <<= logRad;
    }
    scale_out = dnscl-upscl;
  }
};

/* class "ac_fft_dif_r2pX_bfp_inpl"
 *
 * class has member functions--
 * void run( )    -- Implements FFT C-View
 * void coverAssert( )              -- Cover & Assert core contains basic Asserts and cover Conditions
 * ac_int shiftcr( )                -- Bitwise Circular Shift
 * ac_int bitrevint( )              -- Bit-Rversal ac_int
 */

template < int N_FFT, int RADIX, int ORDER, int TWIDDLE_PREC, int DIF_D_P, int DIF_D_I, class out_struct > 
class ac_fft_dif_r2pX_bfp_inpl
{
private:
  typedef ac_fixed < DIF_D_P, DIF_D_I, true > dif_fxp_data;
  typedef ac_complex < dif_fxp_data > dif_input, dif_output;

  ac_fft_dif_r2pX_bfp_inpl_core < N_FFT, RADIX, ORDER, TWIDDLE_PREC, DIF_D_P, DIF_D_I > fft;
  dif_input bank[RADIX][N_FFT / RADIX];

  /*
   * coverAssert() used for basic template assert condition. This helps to validate if object
   * of this class created in user code has right set of parameters defined for it. Code will assert
   * during compile time if incorrect template values are used.
   */

  void coverAssert () {
#ifdef ASSERT_ON
    AC_ASSERT(((N_FFT == 2) || (N_FFT == 4) || (N_FFT == 8) || (N_FFT == 16) || (N_FFT == 32) || (N_FFT == 64) || (N_FFT == 128) ||
               (N_FFT == 256) || (N_FFT == 512) || (N_FFT == 1024) || (N_FFT == 2048) || (N_FFT == 4096)), "N_FFT is not a power of two");
    AC_ASSERT(TWIDDLE_PREC <= 32, "Twiddle bitwidth greater than 32");
    AC_ASSERT(TWIDDLE_PREC >= 2,  "Twiddle bitwidth lesser than 2");
    AC_ASSERT(DIF_D_P >= DIF_D_I,  "Stage integer width lesser than bitwidth");
#endif
#ifdef COVER_ON
    cover (TWIDDLE_PREC <= 5);
#endif
  }

  template < int n, int shft > ac_int < n, 0 > shiftcr (ac_int < n, 0 > &input) {
    const int logn = ac::log2_ceil < n >::val;
    ac_int < n, 0 > output;
    ac_int < logn + 1, 0 > ckr;
#pragma unroll yes
    CIRCULAR_SHIFT: for (ac_int < logn + 1, 0 > itr = 0;; itr++) {
      ckr = (itr + shft) % n;
      output[itr] = input[ckr];

      if (itr == n - 1) { break; }
    }
    return output;
  }

  template < int a > ac_int < a, 0 > bitrevint (ac_int < a, 0 > &Num) {
    ac_int < a, 0 > Num_br;

#pragma unroll yes
    BITREVERSE: for (ac_int < a + 1, 0 > itrator = 0; itrator < a; itrator++) {
      Num_br[a - 1 - itrator] = Num[itrator]; // bitreversal ac_int
    }
    return Num_br;
  }
public:

  ac_fft_dif_r2pX_bfp_inpl () { // Constructor
    ac::init_array < AC_VAL_DC > (&bank[0][0], N_FFT);
  }

  /*
   * Function 'run()' is top of generic RADIX  FFT (DIF) In-place architecture. This  can be called in  user  code
   * by instantiating the object of the class with specific template parameters. It uses ac_channels to explicitly stream
   * design's inputs and outputs in catapult synthesis.
   */

#pragma hls_design interface
  void run(ac_channel < dif_input > &inst, ac_channel < out_struct > &outst) {
    void coverAssert ();

    const int logN = ac::log2_ceil < N_FFT >::val;
    const int logRad = ac::log2_ceil < RADIX >::val;
    const int mixR = logN % logRad;
    const ac_int < logN - logRad + 1, 0 > fix_zero_wire_wdth = ((1 << (logN - logRad)) - 1);
    const int RADIX_R = (logRad - mixR) % logRad;
    const ac_fixed < DIF_D_P, DIF_D_I, true > VAL_MAX_sft1 = (~((ac_fixed < DIF_D_P, DIF_D_I, true >) (1 << (DIF_D_I - 1))));

    ac_int < ac::log2_ceil < DIF_D_P >::val, 1 > upscl = 0;
    dif_input read_temp;
    dif_fxp_data dynamic_range = 0;
    out_struct outp_temp;

    ac_int < ac::log2_ceil < ac::log2_ceil < N_FFT >::val >::val + 1, 1 > scale_out_temp;
#ifndef __SYNTHESIS__
    while (inst.available (N_FFT))
#endif
    {
#pragma hls_pipeline_init_interval 1
      INPUT_LOOP: for (int k = 0; k < RADIX; k++) {
        INPUT_LOOP_ADDRESSING: for (int j = 0; j < N_FFT / RADIX; j++) {
          ac_int < logRad, 0 > bank_temp, bank_add;
          bank_temp = k;
          bank_add = shiftcr < logRad, RADIX_R > (bank_temp);
          read_temp = inst.read ();
          dynamic_range = dynamic_range | (read_temp.r () >= 0 ? read_temp.r () : (~(read_temp.r ()))) | (read_temp.i () >= 0 ? read_temp.i () : (~(read_temp.i ())));
          bank[bank_add][j] = read_temp;
        }
      }

      int f_onebit;

#pragma unroll yes
      FIRST_ONE_LOOP: for (f_onebit = DIF_D_P - 1; f_onebit > 0; f_onebit--)
        if (dynamic_range[f_onebit] == 1) { break; }

      upscl = DIF_D_P - f_onebit - logRad - 2;

      if (upscl < 0) { upscl = 0; }

      fft.fftDifR2pXInplCore (bank, upscl, scale_out_temp);
      outp_temp.expo = scale_out_temp;
      ac_int < logN, 0 > out_add = 0;
      ac_int < logN - logRad + 1, 0 > mem_intr, mem_add;
      ac_int < logRad, 0 > m_no_msb = 0, bank_add;
      if (ORDER == 1) {
        /*
         *Loop Generate Address to Fetch Output in Natural Order
         */
#pragma hls_pipeline_init_interval 1
        OUTPUT_LOOP_NATURAL: for (int m = 0; m < RADIX; m++) {
          m_no_msb = m;
          bank_add = bitrevint < logRad > (m_no_msb);
          m_no_msb = bank_add;
          OUTPUT_LOOP_NATURAL_ADDRESS: for (int i = 0; i < N_FFT / RADIX; i++) {
#pragma unroll yes
            OUTPUT_LOOP_NATURAL_SLICE: for (int p = 0; p <= (logN - logRad); p = (p + logRad)) {
              out_add.set_slc (p, m_no_msb);
            }
            mem_intr = i;
            mem_intr = bitrevint < logN - logRad + 1 > (mem_intr);
            mem_intr = mem_intr >> 1;
            mem_add = ((fix_zero_wire_wdth) & (mem_intr ^ out_add));
            outp_temp.manti = bank[bank_add][mem_add];
            outst.write (outp_temp);
          }
        }
      } else {
        /*
         *Loop Generate Address to Fetch Output in Natural Order
         */
#pragma hls_pipeline_init_interval 1
        OUTPUT_LOOP_BITREVERSE: for (int i = 0; i < N_FFT / RADIX; i++) {
          OUTPUT_LOOP_BITREVERSE_ADDRESS: for (int m = 0; m < RADIX; m++) {
            m_no_msb = m;
#pragma unroll yes
            OUTPUT_LOOP_BITREVERSE_SLICE: for (int p = 0; p <= (logN - logRad); p = (p + logRad)) {
              out_add.set_slc (p, m_no_msb);
            }
            mem_add = ((fix_zero_wire_wdth) & (i ^ out_add));
            bank_add = m;
            outp_temp.manti = bank[bank_add][mem_add];
            outst.write (outp_temp);
          }
        }
      }
#ifndef __SYNTHESIS__
      break;
#endif
    }
  }
};

#endif

