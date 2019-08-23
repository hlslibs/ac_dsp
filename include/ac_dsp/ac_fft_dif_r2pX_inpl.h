/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) DSP Library                                        *
 *                                                                        *
 *  Software Version: 3.2                                                 *
 *                                                                        *
 *  Release Date    : Fri Aug 23 10:38:50 PDT 2019                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.2.0                                               *
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
// File: ac_fft_dif_r2pX_inpl.h
//
// Description:
//  The ac_fft_dif_r2pX_inpl class serves as a C++ interface to the core class. The
//  ac_fft_dif_r2pX_inpl_core class is generic and intended to work well even with SystemC
//  implementations in addition to the supplied C++ implementation. The core class
//  instantiates the butterfly as an object of the 'ac_fft_dif_r2pX_inpl_butterfly' class,
//  handles computation of the FFT Flow Graph and also fetches and writes data.
//  The butterfly class implements the Radix-engine compromise, and the output of it is
//  hard-coded to scale-down by half.
//
// Nomenclature:
//       ac_fft_dif_r2pX_inpl
//               |   |   |
//               |   |   -------- In-place Architecture
//               |   ------------ Radix-2^x
//               ---------------- Decimation in Frequency
//
// Order of Input/Output:
//  Input  -- Natural
//  Output -- ORDER = 0 --> Bit reversed output
//                  = 1 --> Natural output
//
// Usage:
//  A sample testbench and its implementation look like this:
//  #include <ac_fft_dif_r2pX_inpl.h>
//
//  #include <mc_scverify.h>
//
//  CCS_MAIN(int arg, char **argc)
//  {
//    // Initialize object of FFT class with 4 FFT points, Radix = 2, bit-reversed output,
//    // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2.
//    ac_fft_dif_r2pX_inpl < 4, 2, 0, 19, 18, 2 > fft_design1;
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

#ifndef _INCLUDED_AC_FFT_DIF_R2PX_INPL_H_
#define _INCLUDED_AC_FFT_DIF_R2PX_INPL_H_

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
// Class: ac_fft_dif_r2pX_inpl_butterfly
// Description: Templatized class for a butterfly
//-------------------------------------------------------------------------

template <unsigned N_FFT, int RADX, class fix_p, class com_p, class complex_round, class com_rnd_ext, class com_mult, class com_tw >
class ac_fft_dif_r2pX_inpl_butterfly
{
public:
  //-------------------------------------------------------------------------
  // Member Function: compute
  // Description: Computes butterfly function according to given parameters
  // x: Input/Output vector.
  // w: Twiddle vector.
  // scale_fac: Scaling Factor
  //
  void compute (int stage_n, com_p x[RADX], const com_tw w[RADX], const int scale_fac) {

    com_rnd_ext xt[RADX];

    complex_round tmp_out;
    com_tw tw[RADX];

#pragma unroll yes
    TWIDDLE_ROUNDING_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      xt[rad_itr] = x[rad_itr];
      tw[rad_itr] = (com_tw) w[rad_itr];
    }

    fft_dif_r2pX_inpl_radix_butterfly (stage_n, xt, scale_fac);
    multRescale (xt, tw, scale_fac);

#pragma unroll yes
    OUTPUT_ASSINGMENT_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      x[rad_itr] = xt[rad_itr];
    }
  }

private:
  //-------------------------------------------------------------------------
  // Member Function: fft_dif_r2pX_inpl_radix_butterfly
  // Description: Implements Radix Engine, and carries out computation as
  // Radix-2 flow graph.
  //
  void fft_dif_r2pX_inpl_radix_butterfly (int stage_n, com_rnd_ext x[RADX], const int scale_fac) {
    const int logN = ac::log2_ceil < N_FFT >::val;
    const int logRad = ac::log2_ceil < RADX >::val;
    const int mixR = logN % logRad;
    com_tw Rw;

#include "twiddlesR_64bits.h"

#pragma unroll yes
    RADIX_STAGE_LOOP: for (ac_int < logRad, 0 > Rstage = logRad - 1; Rstage >= 0; Rstage--) {
      com_rnd_ext a[RADX];

#pragma unroll yes
      INPUT_REG_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
        a[rad_itr] = x[rad_itr];
      }

#pragma unroll yes
      RADIX_BUTTERFLY_LOOP: for (ac_int < logRad, 0 > B_count = 0; B_count < RADX / 2; B_count++) {
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
          temp_tw.r() = 1;
          temp_tw.i() = 0;
        }
        x[addrs2] = (complex_round) (temp_x * temp_tw);
        x[addrs1] = (complex_round) (x[addrs1]);
      }

#pragma unroll yes
      RESCALE_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
        com_rnd_ext temp;

        bool ch_scale = 1 & (scale_fac >> (logRad - 1 - Rstage));

        if (ch_scale && (Rstage != 0)) {
          temp.r() = (x[rad_itr].r() >> 1);
          temp.i() = (x[rad_itr].i() >> 1);
        } else {
          temp = x[rad_itr];
        }

        x[rad_itr] = temp;
      }

      if (Rstage == 0)
      { break; }
    }
  }

  //-------------------------------------------------------------------------
  // Member Function: multRescale
  // Description: Computes multiplication of twiddle factor and rescales
  // output of multiplication.
  //
  void multRescale (com_rnd_ext x[RADX], com_tw y[RADX], const int scale_fac) {
    const int logRad = ac::log2_ceil < RADX >::val;
    com_mult a[RADX], b[RADX], tx[RADX];
    complex_round temp;
    bool ch_scale = 1 & (scale_fac >> (logRad - 1));;

#pragma unroll yes
    TWIDDLE_MULT_LOOP: for (int rad_itr = 1; rad_itr < RADX; rad_itr++) {
      a[rad_itr] = x[rad_itr];
      b[rad_itr] = y[rad_itr];
      tx[rad_itr] = a[rad_itr] * b[rad_itr];
      if (ch_scale) {
        temp.r() = (tx[rad_itr].r() >> 1);
        temp.i() = (tx[rad_itr].i() >> 1);
      } else {
        temp = tx[rad_itr];
      }
      x[rad_itr] = temp;
    }

    temp = x[0];

    if (ch_scale) {
      x[0].r() = (temp.r() >> 1);
      x[0].i() = (temp.i() >> 1);
    } else {
      x[0] = temp;
    }
  }
};

//=========================================================================
// Class: ac_fft_dif_r2pX_inpl_core
// Description: Core class for FFT computation.
//-------------------------------------------------------------------------

template < unsigned N_FFT, int RADIX, int ORDER, int TWID_PREC, int DIF_D_P, int DIF_D_I >
class ac_fft_dif_r2pX_inpl_core
{
private:
  // Typedefs for public function args declared first, to avoid compile-time errors.
  typedef ac_fixed < DIF_D_P, DIF_D_I, true > dType;
  typedef ac_complex < dType > cx_dType;

public:
  //-------------------------------------------------------------------------
  // Member Function: bitrevint
  // Description: Carries out bit-reversal on ac_int variables.
  //
  template < int a > ac_int < a, 0 > bitrevint (ac_int < a, 0 > &Num) {
    ac_int < a, 0 > Num_br;

#pragma unroll yes
    BITREVERSAL_CORE: for (ac_int < a + 1, 0 > itrator = 0; itrator < a; itrator++) {
      Num_br[a - 1 - itrator] = Num[itrator];
    }
    return Num_br;
  }

  //-------------------------------------------------------------------------
  // Member Function: fftDifR2pXInplCore
  // Description: Implements core functionality of the FFT.
  //
  void fftDifR2pXInplCore (cx_dType bank[RADIX][N_FFT / RADIX]) {

#include "twiddles_64bits.h"

    const int logN_a = ac::log2_ceil < N_FFT >::val;
    const int logRad = ac::log2_ceil < RADIX >::val;
    const int mixR = logN_a % logRad;
    const int RADIX_R = (logRad - mixR) % logRad;
    const int Nstage_a = logN_a / logRad;
    const int modi_sh = (Nstage_a + (mixR != 0)) * logRad;
    const unsigned N_FFT_m = 1 << (modi_sh);
    const int logN = ac::log2_ceil < N_FFT_m >::val;
    const int Nstage = logN / logRad;

    cx_dType data[RADIX];

    ac_int < logN + 1, false > n1, n2;
    ac_int < logN + 1, false > idx;
    ac_int < logN - logRad + 1, 0 > bank_add_gen = 0;
    ac_int < logN, 0 > xadd = 0;
    ac_int < logN - logRad + 1, 0 > id_bank[RADIX], cur_id_bank, fet_id_data;
    const ac_int < logN - logRad + 2, 0 > fix_zero_wire_wdth = ((1 << (logN - logRad)) - 1);

    n1 = N_FFT_m;
    n2 = 0;
    idx = 1;
    cx_tType twd, tw, tw_vac[RADIX];

    // Radix-2^x engine instantiation

    ac_fft_dif_r2pX_inpl_butterfly < N_FFT, RADIX, dType, cx_dType, cx_dround, cx_b_dround, cx_mType, cx_tType > btrfly;

#pragma hls_pipeline_init_interval 1
    STAGE_LOOP: for (int i = Nstage - 1; i >= 0; i--) {
      bank_add_gen = 0;
      n2 = n1;
      n1 >>= logRad;
      SEGMENT_LOOP: for (int j = 0; j < N_FFT / RADIX; j++) {
        int k = j;

        BUTTERFLY_LOOP: for (int kk = 0; kk < N_FFT / RADIX; kk++) {

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
            data[bank_add_fet] = bank[mn][fet_id_data];
          }
          ac_int < (logN + 1), false > n, n_vac[RADIX];
          n = (j * idx);

#pragma unroll yes
          BITREVERCE_TWIDDLE_ADDRESS: for (ac_int < logRad + 1, 0 > n_itr = 0; n_itr < RADIX; n_itr++) {
            ac_int < logRad, 0 > n_nmsb = n_itr;
            n_nmsb = bitrevint(n_nmsb);
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
            tw.r() = ((1 & ((n << 2) >> (logN_a - 1))) | (1 & ((n << 1) >> (logN_a - 1)))) ? (tType) (-twd.r()) : (tType) (twd.r());
            tw.i() = ((!(1 & ((n << 2) >> (logN_a - 1)))) | (1 & ((n << 1) >> (logN_a - 1)))) ? (tType) (twd.i()) : (tType) (-twd.i());
            tw = ((1 & ((n << 2) >> (logN_a - 1))) ^ (1 & ((n << 1) >> (logN_a - 1)))) && ((logN_a - 1) >= 2) ? (cx_tType) (J * tw.conj ()) : (cx_tType) tw;
            tw = n[logN_a - 1] ? (cx_tType) (-tw) : (cx_tType) (tw);
            tw_vac[tw_itr] = tw;
          }

          // Note: this design is hard-coded to scale all outputs by 1/2. If you do not want that to happen, substitute the "N_FFT - 1" in the equation below with an
          // integer whose bit pattern encodes the scaling of all stages, where the LSB is scaling of the output stage and MSB is scaling of the input stage.
          btrfly.compute (i - Nstage + 1, data, tw_vac, ((RADIX - 1) & (((N_FFT - 1) << RADIX_R) >> ((Nstage - 1 - i) * logRad))));

#pragma unroll yes
          PUTTING_DATA_IN_BANKS: for (int amn = 0; amn < RADIX; amn++) {
            cur_id_bank = id_bank[amn] & ((N_FFT / RADIX) - 1);
            int data_add = amn ^ (bank_add_gen.template slc < logRad > (i > 0 ? (i - 1) * logRad : logN));

            bank[amn][cur_id_bank] = data[data_add];
          }
          k += n2;
          bank_add_gen += (1 << ((i) * logRad));

          if ((bank_add_gen >= N_FFT / RADIX)) {
            bank_add_gen = bank_add_gen & ((N_FFT / RADIX) - 1);
            bank_add_gen += 1;
          }

          if (k >= N_FFT - 1)
          { break; }
        }
        if (j == n1 - 1)
        { break; }

      }
      idx <<= logRad;
    }
  }
  
private:
  typedef ac_fixed < TWID_PREC, 2, true, AC_RND_INF > tType;
  typedef ac_complex < tType > cx_tType;
  typedef ac_fixed < DIF_D_P + 1, DIF_D_I + 1, true, AC_RND, AC_SAT > dround;
  typedef ac_complex < dround > cx_dround;
  typedef ac_fixed < DIF_D_P + TWID_PREC - 1, 1 + DIF_D_I, true > mType;
  typedef ac_complex < mType > cx_mType;
  typedef ac_fixed < DIF_D_P + 1, DIF_D_I + 1, true > b_dround;
  typedef ac_complex < b_dround > cx_b_dround;
};

//==================================================================================
// Class: ac_fft_dif_r2pX_inpl
// Description:
// Top-level class of design, instantiated by testbench.
// HLS Interface: run()
//----------------------------------------------------------------------------------

template < unsigned N_FFT, int RADIX, int ORDER, int TWIDDLE_PREC, int DIF_D_P, int DIF_D_I >
class ac_fft_dif_r2pX_inpl
{
private:
  // Typedefs for public function args declared first, to avoid compile-time errors.
  typedef ac_fixed < DIF_D_P, DIF_D_I, true > dif_fxp_data;
  typedef ac_complex < dif_fxp_data > dif_input, dif_output;
  
public:
  //-------------------------------------------------------------------------
  // Constructor
  //
  ac_fft_dif_r2pX_inpl () {
    ac::init_array < AC_VAL_DC > (&bank[0][0], N_FFT);
  }

  //----------------------------------------------------------------------
  // Member Function: run
  // Description: the top-level interface to the FFT.
  //
#pragma hls_design interface
  void CCS_BLOCK(run)(ac_channel < dif_input > &inst, ac_channel < dif_output > &outst) {
    coverAssert();
    const int logN = ac::log2_ceil < N_FFT >::val;
    const int logRad = ac::log2_ceil < RADIX >::val;
    const int mixR = logN % logRad;
    const ac_int < logN - logRad + 1, 0 > fix_zero_wire_wdth = ((1 << (logN - logRad)) - 1);
    const int RADIX_R = (logRad - mixR) % logRad;

#pragma hls_pipeline_init_interval 1
    INPUT_LOOP: for (int k = 0; k < RADIX; k++) {

      INPUT_LOOP_ADDRESSING: for (int j = 0; j < N_FFT / RADIX; j++) {
        ac_int < logRad, 0 > bank_temp, bank_add;
        bank_temp = k;
        bank_add = shiftcr <RADIX_R> (bank_temp);
        bank[bank_add][j] = inst.read ();
      }
    }

    fft.fftDifR2pXInplCore (bank);

    ac_int < logN, 0 > out_add = 0;
    ac_int < logN - logRad + 1, 0 > mem_intr, mem_add;
    ac_int < logRad, 0 > m_no_msb = 0, bank_add;

    if (ORDER == 1) {
      // This loop Generates Addresses to Fetch Output in Natural Order
#pragma hls_pipeline_init_interval 1
      OUTPUT_LOOP_NATURAL: for (int m = 0; m < RADIX; m++) {
        m_no_msb = m;
        bank_add = fft.bitrevint(m_no_msb);
        m_no_msb = bank_add;
        OUTPUT_LOOP_NATURAL_ADDRESS: for (int i = 0; i < N_FFT / RADIX; i++) {
#pragma unroll yes
          OUTPUT_LOOP_NATURAL_SLICE: for (int p = 0; p <= (logN - logRad); p = (p + logRad)) {
            out_add.set_slc (p, m_no_msb);
          }
          mem_intr = i;
          mem_intr = fft.bitrevint(mem_intr);
          mem_intr = mem_intr >> 1;
          mem_add = ((fix_zero_wire_wdth) & (mem_intr ^ out_add));
          outst.write (bank[bank_add][mem_add]);
        }
      }
    } else {
      // This loop Generates Addresses to Fetch Output in Bit-reversed Order
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
          outst.write (bank[bank_add][mem_add]);
        }
      }
    }
  }

  //-------------------------------------------------------------------------
  // Member Function: coverAssert
  // Description: used for basic template assert condition. This helps to 
  // validate if object of this class created by the user has right set of
  // parameters defined for it. Code will assert during compile-time if incorrect
  // template values are used.
  //
  void coverAssert () {
#ifdef ASSERT_ON
    static_assert(N_FFT == 2   || N_FFT == 4   || N_FFT == 8    || N_FFT == 16   || N_FFT == 32   || N_FFT == 64 || N_FFT == 128 ||
                  N_FFT == 256 || N_FFT == 512 || N_FFT == 1024 || N_FFT == 2048 || N_FFT == 4096, "N_FFT is not a power of two");
    static_assert(RADIX == 2 || RADIX == 4 || RADIX == 8 || RADIX == 16 || RADIX == 32 || RADIX == 64 || RADIX == 128, "Radix Value not supported");
    static_assert(TWIDDLE_PREC <= 32, "Twiddle bitwidth greater than 32");
    static_assert(TWIDDLE_PREC >= 2,  "Twiddle bitwidth lesser than 2");
    static_assert(DIF_D_P >= DIF_D_I,  "Stage integer width lesser than bitwidth");
#endif
#ifdef COVER_ON
    cover (TWIDDLE_PREC <= 5);
#endif
  }
  
private:
  ac_fft_dif_r2pX_inpl_core < N_FFT, RADIX, ORDER, TWIDDLE_PREC, DIF_D_P, DIF_D_I > fft;
  dif_input bank[RADIX][N_FFT / RADIX];

  //-------------------------------------------------------------------------
  // Member Function: shiftcr
  // Description: Implements bitwise circular shift.
  //
  template < int shft, int n >
  ac_int < n, 0 > shiftcr (ac_int < n, 0 > &input) {
    const int logn = ac::log2_ceil < n >::val;

    ac_int < n, 0 > output;
    ac_int < logn + 1, 0 > ckr;
#pragma unroll yes
    CIRCULAR_SHIFT: for (ac_int < logn + 1, 0 > itr = 0; itr < n; itr++) {
      ckr = (itr + shft) % n;
      output[itr] = input[ckr];
    }
    return output;
  }
};

#endif

