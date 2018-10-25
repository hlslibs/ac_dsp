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
// File: ac_fft_dif_r2pX_dyn_inpl.h
//
// Description:
//    Nomenclature:
//                           ac_fft_dif_r2pX_dyn_inpl
//                         /    /    |   \     \      \                              
//                        /   FFT    |  R-2^x   \      \                             
//                  C-view           |          Dynamic  In-place
//                        Decimation in Frequency
//
//    Organization:
//     The ac_fft_dif_r2pX_dyn_inpl class serves as a C++ interface to the core class. The
//     ac_fft_dif_r2pX_dyn_inpl_core class is generic and intended to work well even with SystemC
//     implementations in addition to the supplied C++ implementation. The core class
//     instantiates the butterfly as an object of the 'ac_fft_dif_r2pX_dyn_inpl_butterfly' class,
//     handles computation of the FFT Flow Graph and also fetches and writes data.
//     The butterfly class implements the Radix-engine compromise, and the output of it is
//     hard-coded to scale-down by half.
//     Order of Input/Output:
//     Input -- Natural
//     Output-- ORDER = 0 --> Bit reversed output
//                    = 1 --> Natural output
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_fft_dif_r2pX_dyn_inpl.h>
//
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Initialize object of FFT class with 4 FFT points, Radix = 2, bit-reversed output,
//      // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2.
//      ac_fft_dif_r2pX_dyn_inpl < 4, 2, 0, 19, 18, 2 > fft_design1;
//      typedef ac_complex<ac_fixed<18, 2, true, AC_TRN, AC_WRAP> > IO_type;
//      // Initialize channels for input, output, dynamic reduction and number of instances.
//      ac_channel<IO_type> input;
//      ac_channel<IO_type> output;
//      ac_channel<ac_int<ac::log2_ceil<ac::log2_ceil<4>::val-ac::log2_ceil<2>::val+1>::val,0> > dynstream;
//      ac_channel<ac_int<ac::log2_ceil<4>::val-ac::log2_ceil<2>::val,0> > inststream;
//      // Write 4 inputs to the input port, and set dynamic reduction = 0, no of instances = 1.
//      input.write(IO_type( .125, .25));
//      input.write(IO_type( .375, .50));
//      input.write(IO_type(1.525, .75));
//      input.write(IO_type( .125, .75));
//      dynstream.write(0);
//      inststream.write(1);
//      // Call the top-level function.
//      fft_design1.run(input, output, dynstream, inststream);
//
//      CCS_RETURN(0);
//    }
//
// Notes:
//    Attempting to call the function with a type that is not implemented will result
//    in a compile error.
//    Currently, the block only accepts signed ac_complex<ac_fixed> inputs and outputs which use AC_TRN
//    and AC_WRAP as their rounding and overflow modes.
//    Dynamic Behaviour : The FFT Architecture is Dynamic in No of FFT Points and Instances
//    to be Calculated in one call.
//
//*********************************************************************************************************

#ifndef _INCLUDED_AC_DIF_R2PX_DYN_INPL_H_
#define _INCLUDED_AC_DIF_R2PX_DYN_INPL_H_

#include <ac_fixed.h>
#include <ac_complex.h>
#include <ac_channel.h>
#ifndef __SYNTHESIS__
#include <iostream>
using namespace std;
#endif
const double pi = 3.141592653589793;

//************************************************************************************
// Description:
// class "ac_fft_dif_r2pX_dyn_inpl_butterfly"
//
// class has member functions--
// fftDifR2pXDynInplRadixButterfly( ) -- Radix Engine in Radix-2 Butterfly Structure
// multResacle( ) -- Multiply with twiddle factor and scale by 1/2
// compute( )     -- Compute FFT butterfly call
//
//************************************************************************************

template < unsigned N_FFT, int RADX, class fix_p, class com_p, class complex_round, class com_rnd_ext, class com_mult, class com_tw > 
class ac_fft_dif_r2pX_dyn_inpl_butterfly
{
  // Declaration of Member Functions
private:
  // function 'fftDifR2pXDynInplRadixButterfly' performs Dynamic Radix Engine computation as Radix-2 flow graph
  void fftDifR2pXDynInplRadixButterfly (ac_int < ac::log2_ceil < RADX >::val, 0 > RADIX_R_dyn, int stage_n, com_rnd_ext x[RADX], com_rnd_ext y[RADX], const int scale_fac) {
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
      INPUT_REG_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
        a[rad_itr] = x[rad_itr];
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
        RADIX_ADD_GEN_LOOP:for (int bit_count = 0; bit_count < logRad; bit_count++) {
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

        if ((((stage_n != 0)) || (!(Rstage >= (logRad - RADIX_R_dyn))) || (RADIX_R_dyn == 0))) {
          x[addrs1] = a[addrs1] + a[addrs2];
          x[addrs2] = a[addrs1] - a[addrs2];
        } else {
          x[addrs1] = a[addrs1];
          x[addrs2] = a[addrs2];
        }
        com_mult temp_x, temp_tw;

        temp_x = x[addrs2];

        if ((((stage_n != 0)) || (!(Rstage >= (logRad - RADIX_R_dyn))) || (RADIX_R_dyn == 0))) {
          temp_tw = Rw;
        } else {
          temp_tw.r () = 1;
          temp_tw.i () = 0;
        }
        x[addrs2] = (com_rnd_ext) (temp_x * temp_tw);
        x[addrs1] = (com_rnd_ext) (x[addrs1]);
      }

#pragma unroll yes
      RESCALE_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
        com_rnd_ext temp;
        bool ch_scale = 1 & (scale_fac >> (logRad - 1 - Rstage));
        if (ch_scale && (Rstage != 0)) {
          temp.r () = (x[rad_itr].r () >> 1);
          temp.i () = (x[rad_itr].i () >> 1);
        } else {
          temp = x[rad_itr];
        }
        x[rad_itr] = temp;
      }

      if (Rstage == 0) { break; }
    }

#pragma unroll yes
    OUTPUT_REG_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      y[rad_itr] = x[rad_itr];
    }
  }

  // function 'multRescale' computes multiplication of twiddle factor and Rescales the Output of Multiplication
  // scale_fac = 0 No scaling
  // scale_fac = 1 Scale down by 1/2
  void multRescale (com_rnd_ext x[RADX],com_rnd_ext yo[RADX], com_tw y[RADX], const int scale_fac) {
    const int logRad = ac::log2_ceil < RADX >::val;
    com_mult a[RADX], b[RADX], tx[RADX];
    complex_round temp;
    com_rnd_ext temps;
    bool ch_scale = 1 & (scale_fac >> (logRad - 1));

#pragma unroll yes
    TWIDDLW_MULT_LOOP: for (int rad_itr = 1; rad_itr < RADX; rad_itr++) {
      a[rad_itr] = x[rad_itr];
      b[rad_itr] = y[rad_itr];
      tx[rad_itr] = a[rad_itr] * b[rad_itr];
      if (ch_scale) {
        temp.r () = (tx[rad_itr].r () >> 1);
        temp.i () = (tx[rad_itr].i () >> 1);
      } else {
        temp = tx[rad_itr];
      }
      yo[rad_itr] = temp;
    }

    temps = x[0];
    if (ch_scale) {
      temp.r () = (temps.r () >> 1);
      temp.i () = (temps.i () >> 1);
    } else {
      temp = temps;
    }
    yo[0] = temp;
  }

public:
  void compute (ac_int < ac::log2_ceil < RADX >::val, 0 > RADIX_R_dyn, int stage_n, com_p x[RADX], com_p y[RADX], const com_tw w[RADX], const int scale_fac) {

    com_rnd_ext xt[RADX];
    com_rnd_ext xyt[RADX];
    com_rnd_ext yt[RADX];

    complex_round tmp_out;
    com_tw tw[RADX];

#pragma unroll yes
    TWIDDLE_ROUNDING_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      xt[rad_itr] = x[rad_itr];
      tw[rad_itr] = (com_tw) w[rad_itr];
    }

    fftDifR2pXDynInplRadixButterfly (RADIX_R_dyn, stage_n, xt, xyt, scale_fac);

    multRescale (xyt,yt, tw, scale_fac);

#pragma unroll yes
    OUTPUT_ASSINGMENT_LOOP: for (int rad_itr = 0; rad_itr < RADX; rad_itr++) {
      y[rad_itr] = yt[rad_itr];
    }
  }

};

//***********************************************************************************
// Description:
//
//  'ac_fft_dif_r2pX_dyn_inpl_core' class has member functions--
//  fftDifR2pXDynInplCore() -- Implements core functionality the FFT.
//  bitrevint() -- Bitrevarsal ac_int
//
//***********************************************************************************

template < unsigned N_FFT, int RADIX, int ORDER, int TWID_PREC, int DIF_D_P, int DIF_D_I >
class ac_fft_dif_r2pX_dyn_inpl_core
{
public:
  // Type definitions for Multipliers,Accumulator and stage variable for fft
  typedef ac_fixed < TWID_PREC, 2, true, AC_RND_INF > tType;
  typedef ac_complex < tType > cx_tType;
  typedef ac_fixed < DIF_D_P, DIF_D_I, true > dType;
  typedef ac_complex < dType > cx_dType;
  typedef ac_fixed < DIF_D_P, DIF_D_I, true, AC_RND, AC_SAT > dround;
  typedef ac_complex < dround > cx_dround;
  typedef ac_fixed < DIF_D_P + TWID_PREC - 1, 1 + DIF_D_I, true > mType;
  typedef ac_complex < mType > cx_mType;
  typedef ac_fixed < DIF_D_P + 1, DIF_D_I + 1, true > b_dround;
  typedef ac_complex < b_dround > cx_b_dround;

  template < int a > ac_int < a, 0 > bitrevint (ac_int < a, 0 > &Num) {
    ac_int < a, 0 > Num_br;

#pragma unroll yes
    BITREVERSAL_CORE: for (ac_int < a + 1, 0 > itrator = 0; itrator < a; itrator++) {
      Num_br[a - 1 - itrator] = Num[itrator]; // bitreversal ac_int
    }
    return Num_br;
  }

  void fftDifR2pXDynInplCore (ac_int < ac::log2_ceil < N_FFT >::val + 1, 0 > N_FFT_dyn,
                              ac_int < ac::log2_ceil < N_FFT >::val, 0 > logN_dyn,
                              cx_dType bank[RADIX][N_FFT / RADIX],
                              ac_int < ac::log2_ceil < N_FFT >::val, 0 > Nstage_dyn,
                              ac_int < ac::log2_ceil < ac::log2_ceil < N_FFT >::val - ac::log2_ceil < RADIX >::val + 1 >::val, 0 > &dyn_red,
                              ac_int < ac::log2_ceil < RADIX >::val, 0 > &RADIX_R_dyn, ac_int < ac::log2_ceil < N_FFT >::val - ac::log2_ceil < RADIX >::val, 0 > &instance_count) {
    // FFT computation will be done in 'fftDifR2pXDynInplCore'
#include <twiddles_64bits.h>
    // Actual FFT Constants
    const int logN_a = ac::log2_ceil < N_FFT >::val;
    const int loglogN_a = ac::log2_ceil < logN_a >::val;
    const int logRad = ac::log2_ceil < RADIX >::val;
    const int mixR = logN_a % logRad;
    const int RADIX_R = (logRad - mixR) % logRad;

    // Modified FFT Constants
    const int logN = logN_a + RADIX_R;
    const int Nstage = (logN_a / logRad) + 1;

    cx_dType data_in[RADIX],data_out[RADIX];

    ac_int < logN + 1, false > n1, n2;
    ac_int < logN + 1, false > idx;
    ac_int < logN - logRad + 1, 0 > bank_add_gen = 0;
    ac_int < logN, 0 > xadd = 0;        //
    ac_int < logN - logRad + 1, 0 > id_bank[RADIX], cur_id_bank, fet_id_data;
    ac_int < logN - logRad + 2, 0 > fix_zero_wire_wdth;

    cx_tType twd, tw, tw_vac[RADIX];

    // Radix -X engine instantiation
    ac_fft_dif_r2pX_dyn_inpl_butterfly < N_FFT, RADIX, dType, cx_dType, cx_dround, cx_b_dround, cx_mType, cx_tType > btrfly;
    idx = 1;
    fix_zero_wire_wdth = (((1 << (logN - logRad)) - 1) >> dyn_red);

    // Loop Structure Merged with addressing
    int i = Nstage_dyn;
    int shifter_num = 0;
    int box_mask = ((1 << dyn_red) - 1);

    ac_int < logN_a, 0 > logRxStage = logN_dyn + RADIX_R_dyn;

    /* Printing on std Output */
#ifndef __SYNTHESIS__
    cout << "FFT Max = " << N_FFT << " Radix = " << RADIX << " Dynamic Reduction = " << dyn_red << " N_FFT dynamic = " << N_FFT_dyn << " Instances = " << instance_count << endl;
#endif

    // Loop Info:FFT_STAGE_LOOP
    //
    // Loop executes Stage wise FFT Computation.
    // loop iteration pattern shown below for example  N_FFT=32 and RADIX=4 (Dynamic FFT Points =16 and Instance count = 2 )
    //
    //  ####### Note: Pipelining @ II = 1 #######
    //  #                                       #
    //  #  The loop Can be Pipelined at II = 1  #
    //  #  but in that case user has to give    #
    //  #  specific "ignore mem precedences".   #
    //  #  If generic "ignore mem precedences"  #
    //  #  are given its safe to lose this      #
    //  #  constraint.                          #
    //  #                                       #
    //  #########################################
    //
    // FFT_STAGE_LOOP: for (Compute Stage Calculation)
    //                     {
    //                                                 loop iteration 1                               loop iteration 2
    //                     Initial state of Banks                         Banks after Stage 1                            Banks after Stage 2
    //                     B1      B2      B3       B4                    B1      B2      B3       B4                    B1      B2      B3       B4
    //
    //                     [0]     [4]     [8]     [12]                   [0]     [4]     [8]     [12]                   [0]     [5]     [10]    [15]
    //                     [1]     [5]     [9]     [13]                   [5]     [1]     [13]    [9]                    [4]     [1]     [14]    [11]
    //                     [2]     [6]     [10]    [14]    ----------\    [10]    [14]    [2]     [6]     ----------\    [8]     [13]    [2]     [7]
    //                     [3]     [7]     [11]    [15]    STAGE 1    \   [15]    [11]    [7]     [3]     STAGE 2    \   [12]    [9]     [6]     [3]
    //                     [0]     [4]     [8]     [12]    COMPUTATION/   [0]     [4]     [8]     [12]    COMPUTATION/   [0]     [5]     [10]    [15]
    //                     [1]     [5]     [9]     [13]    ----------/    [5]     [1]     [13]    [9]     ----------/    [4]     [1]     [14]    [11]
    //                     [2]     [6]     [10]    [14]                   [10]    [14]    [2]     [6]                    [8]     [13]    [2]     [7]
    //                     [3]     [7]     [11]    [15]                   [15]    [11]    [7]     [3]                    [12]    [9]     [6]     [3]
    //
    //                     }

    FFT_STAGE_LOOP: for (int iterator = 0; iterator < Nstage; iterator++) {
      i--;
      logRxStage -= logRad;
      bank_add_gen = 0;

      ac_int < logN_a - logRad, 0 > j = 0;
      ac_int < logN_a - logRad, 0 > box = 0;
      ac_int < logN_a - logRad, 0 > uni_counter = 0;

      // Loop Info:BUTTERFLY_CALL_LOOP
      //
      // Loop executes butterfly Computation for one stage.
      // All Memory Addressing is done in this loop.
      // loop iteration pattern shown below for example  N_FFT=32 and RADIX=4 (Dynamic FFT Points =16 and Instance count = 2 )
      //
      // BUTTERFLY_CALL_LOOP: for (Compute Butterfly Calculation)
      //                     {
      //
      //                     Initial state of Banks                      First stage 3rd Butterfly call
      //                     B1      B2      B3       B4     ---------\
      //                                                       FEED    \       D  C  B  A
      //          loop itr 1 [0]     [4]     [8]     [12]     13 9 5 1 /         \ \/ /
      //          loop itr 3 [X]     [X]     [X]     [X]     ---------/           \/\/       <---   Radix-4
      //              .      [2]     [6]     [10]    [14]                         /\/\          Butterfly Engine
      //              .      [3]     [7]     [11]    [15]                        / /\ \
      //          loop itr 2 [0]     [4]     [8]     [12]   /---------          * *  | |
      //          loop itr 4 [1]     [5]     [9]     [13]  /  OUT               \/    \/
      //              .      [2]     [6]     [10]    [14]  \ 5 1 13 9           /\    /\
      //              .      [3]     [7]     [11]    [15]   \---------         D  C  B  A
      //
      //                     }

#pragma hls_pipeline_init_interval 1
      BUTTERFLY_CALL_LOOP: for (ac_int < logN_a - logRad, 0 > uni_counter1 = 0; uni_counter1 < (1 << (logN_a - logRad)); uni_counter1++) {
        ac_int < logN_a - logRad, 0 > uni_counter_dyn_shifted;
        uni_counter_dyn_shifted = (uni_counter >> dyn_red);
        j = uni_counter_dyn_shifted >> shifter_num;
        box = uni_counter & box_mask;

        //  Loop Info:BANK_ADD_GENERATOR
        //
        //  loop calculates memory address for single Butterfly call
        //  loop will kept fully unrolled to calculate required address
        //  in one clock cycle
        //
        //  Address of Coefficients depends on stage and Butterfly No.
        //  example for 3rd Butterfly of N_FFT=32 and RADIX=4 (Dynamic FFT Points =16 and Instance count = 2 )
        //  given below
        //
        // BANK_ADD_GENERATOR: for (Calculate Address)
        //                      {
        //
        //                      Initial state of Banks
        //                      B1      B2      B3       B4     ADDRESS     ---------\
        //                                                     id_bank[]      FEED    \
        //           loop itr 1 [0]     [4]     [8]     [12]    { B1-1       13 9 5 1 /
        //           loop itr 3 [X]     [X]     [X]     [X]       B2-1      ---------/
        //               .      [2]     [6]     [10]    [14]      B3-1
        //               .      [3]     [7]     [11]    [15]      B4-1}
        //           loop itr 2 [0]     [4]     [8]     [12]               /---------
        //           loop itr 4 [1]     [5]     [9]     [13]              /  OUT
        //               .      [2]     [6]     [10]    [14]              \ 5 1 13 9
        //               .      [3]     [7]     [11]    [15]               \---------
        //
        //                      }

#pragma unroll yes
        BANK_ADD_GENERATOR: for (ac_int < logRad + 1, 0 > m = 0; m < RADIX; m++) {
          ac_int < logRad, 0 > m_no_msb = 0;
          xadd = 0;
          m_no_msb = m;
#pragma unroll yes
          BANK_ADD_GEN_NEST: for (int slc_j = logN - logRad; slc_j >= 0; slc_j = (slc_j - logRad)) {
            if (slc_j < logRxStage) { break; }
            xadd.set_slc (slc_j, m_no_msb);
          }
          id_bank[m] = fix_zero_wire_wdth & ((bank_add_gen) ^ xadd);
        }

        //  Loop Info:DATA_FETCH_FROM_BANKS
        //
        //  Loop fetches data for Butterfly Computation
        //  Always Fully Unroll this loop
        //
        // DATA_FETCH_FROM_BANKS: for (FETCH DATA)
        //                      {
        //
        //                     data_in[x] = memory_banks[x]
        //                     data_in[] will be pass to butterfly call
        //                      }

#pragma unroll yes
        DATA_FETCH_FROM_BANKS: for (int mn = 0; mn < RADIX; mn++) {
          int bank_add_fet = mn ^ (bank_add_gen.template slc < logRad > (logRxStage));
          fet_id_data = id_bank[mn] & (((N_FFT / RADIX) - 1) >> dyn_red);
          fet_id_data = (fet_id_data | (box << (logN_a - logRad - dyn_red))) & ((N_FFT / RADIX) - 1);
          data_in[bank_add_fet] = bank[mn][fet_id_data];
        }
        ac_int < (logN + 1), false > n, n_vac[RADIX];

        // if Scope is for these calculations which are the same for one stage
        if (box == 0) {
          n = (j * idx) << dyn_red;
          //  Loop Info:BITREVERSE_TWIDDLE_ADDRESS
          //
          //  Bit-reverse twiddles are multiplied each after Radix-Engine Computation
          //  this loop generates Bitreverse twiddle address
          // BITREVERSE_TWIDDLE_ADDRESS : for (Compute twiddle ROM Address)
          //                      {
          //
          //                       Radix-2             Radix-4                                Radix-8
          //      n_bit_rev[]      0     1        0    2    1    3             0    4    2    6    1    5    3    7
          //
          //                      Twiddle Address = n_bi_rev[]  X  n(multiplication factor )
          //
          //                      }

#pragma unroll yes
          BITREVERSE_TWIDDLE_ADDRESS: for (ac_int < logRad + 1, 0 > n_itr = 0; n_itr < RADIX; n_itr++) {
            ac_int < logRad, 0 > n_nmsb = n_itr;
            n_nmsb = bitrevint < logRad > (n_nmsb);
            if ((i == Nstage_dyn - 1) && (RADIX_R_dyn != 0)) {
              n_vac[n_nmsb] = ((n_itr >> RADIX_R_dyn) * (n + ((n_nmsb >> (logRad - RADIX_R_dyn)) * ((N_FFT / RADIX)))));
            } else {
              n_vac[n_nmsb] = ((n_itr * n) >> RADIX_R_dyn);
            }

          }
          cx_tType J = comx_twiddle (0, 1);

          tw_vac[0] = comx_twiddle (1, 0);
          ac_int < logN, 0 > logNminone = logN_a - 1;

          //  Loop Info:TWIDDLE_VAC_GEN
          //
          //  Loop Fetch twiddle value from ROM
          //  #########   NOTE  ROM  ##########
          //  #                               #
          //  #  Twiddles are mapped to ROM   #
          //  #  only 1/8 Cycle of complex    #
          //  #  exponential is Stored        #
          //  #                               #
          //  #################################
          //  TWIDDLE_VAC_GEN : for (Compute Stage Calculation)
          //                      {
          //                        twiddle values[ ]    =  ROM[ twiddle Address ]
          //                      }
          //

#pragma unroll yes
          TWIDDLE_VAC_GEN:for (ac_int < logRad + 1, 0 > tw_itr = 1; tw_itr < RADIX; tw_itr++) {
            ac_int < logN, false > t;
            n = n_vac[tw_itr] & (N_FFT - 1);
            t = (1 & ((n << 2) >> (logNminone))) ? (ac_int < (logN) + 1, false >) ((((1 << (logNminone)) >> 2)) - (n & ((((1 << (logNminone)) >> 2)) - 1))) : (ac_int < (logN) + 1, false >) (n & ((((1 << (logNminone)) >> 2)) - 1));
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
            tw.r () = ((1 & ((n << 2) >> (logNminone))) | (1 & ((n << 1) >> (logNminone)))) ? (tType) (-twd.r ()) : (tType) (twd.r ());
            tw.i () = ((!(1 & ((n << 2) >> (logNminone)))) | (1 & ((n << 1) >> (logNminone)))) ? (tType) (twd.i ()) : (tType) (-twd.i ());
            tw = ((1 & ((n << 2) >> (logNminone))) ^ (1 & ((n << 1) >> (logNminone)))) && ((logNminone) >= 2) ? (cx_tType) (J * tw.conj ()) : (cx_tType) tw;
            tw = n[logNminone] ? (cx_tType) (-tw) : (cx_tType) (tw);
            tw_vac[tw_itr] = tw;
          }
        }                   // BOX == 0 Loop
        int temp_count = 0;

        // Note: this design is hard-coded to scale all outputs by 1/2. If you do not want that to happen, substitute the "N_FFT - 1" in the equation below with an
        // integer whose bit pattern encodes the desired scaling of all stages, where the LSB is scaling of the output stage and MSB is scaling of the input stage.
        btrfly.compute (RADIX_R_dyn, i - Nstage_dyn + 1, data_in, data_out, tw_vac, ((RADIX - 1) & (((N_FFT - 1) << RADIX_R_dyn) >> ((Nstage_dyn - 1 - i) * logRad))));

        //  Loop Info:PUTTING_DATA_IN_BANKS
        //
        //  Loop executes to put data back to Banks after Butterfly Computation
        //  Always Fully Unrolled loop
        //
        // PUTTING_DATA_IN_BANKS: for (PUT DATA)
        //                      {
        //
        //                     Butterfly_out is Data after Butterfly computation and twiddle multiplication
        //                      memory_banks[x] = Butterfly_out[x]
        //
        //                      }
#pragma unroll yes
        PUTTING_DATA_IN_BANKS: for (int amn = 0; amn < RADIX; amn++) {
          cur_id_bank = id_bank[amn] & (((N_FFT / RADIX) - 1) >> dyn_red);
          cur_id_bank = (cur_id_bank | box << (logN_a - logRad - dyn_red)) & ((N_FFT / RADIX) - 1);
          int data_add = amn ^ (bank_add_gen.template slc < logRad > (i > 0 ? ((int) (logRxStage) - logRad) : logN));
          bank[amn][cur_id_bank] = data_out[data_add];
        }

        if (box == box_mask || ((instance_count != 0) && (box == (instance_count - 1)))) {
          bank_add_gen += (1 << logRxStage);
          if ((bank_add_gen >= ((N_FFT / RADIX) >> dyn_red))) {
            bank_add_gen = bank_add_gen & (((N_FFT / RADIX) - 1) >> dyn_red);
            bank_add_gen += 1;
          }
        }

        if ((instance_count != 0) && (box == (instance_count - 1))) {
          uni_counter += (box_mask - instance_count + 1);
        }

        if (uni_counter == (((N_FFT / RADIX) - 1))) { break; }
        uni_counter++;
      }                       // Universal Counter Loop End

      idx <<= logRad;
      if ((iterator == 0) && (RADIX_R_dyn != 0)) {
        shifter_num += (logRad - RADIX_R_dyn);
      } else {
        shifter_num += logRad;
      }

      if (i == 0) { break; }              // Breaking Stage Loop
    }                           // Stage Loop End
  }
};

#ifdef ASSERT_ON
#include <ac_assert.h>
#endif
#ifdef COVER_ON
#include <ac_cover.h>
#endif

//**************************************************************************************************************************
//
// class "ac_fft_dif_r2pX_dyn_inpl"
//
// class has member functions--
// void run( )                      -- Implements FFT C-View
// void coverAssert( )              -- Cover & Assert core contains basic Asserts and cover Conditions
// ac_int shiftcr( )                -- Bitwise Circular Shift
// ac_int bitrevint( )              -- Bit-Rversal ac_int
//
//
//  Template Info:
//  template < unsigned N_FFT, int RADIX, int ORDER, int TWIDDLE_PREC, int DIF_D_P, int DIF_D_I >
//
//  N_FFT       : Maximum No. of FFT Points can Calculated in one call
//  ORDER       : Order of Output 0 = Bitreversed, 1 = Natural
//  TWIDDLE_PREC: Twiddle Bit Precision
//  DIF_D_P     : Bit-precision of Calculation
//  DIF_D_I     : Integer Part in Bits
//  Port Info:
//  run (ac_channel inst, ac_channel outst, ac_channel dyn_red_chn, ac_channel ins_count);
//
//  Input Ports:
//  Three Input Ports
//  1. Input Data Port   : Port is used for Data Input of given width. Order of Data is expected to be Sequential.
//  2. Dynamic Reduction : Port is used for Change in FFT Points Dynamically for value = 1 on this Port No. FFT
//                         Points will reduce to N_FFT/2 similarly for =2, N_FFT/4 and so on till Dynamic FFT
//                         Points reaches value equals to RADIX
//  3. No of Dynamic FFT : Port is used for No. of Dynamic FFT Instance. Maximum value of also depends
//        Instances        on value Dynamic Reduction eg. Dynamic Reducton = 2 then, FFT Point = N_FFT/4 it means
//                         RAM is capable of N_FFT( Maximum No. of FFT Points ) only N_FFT/4 is used thus instance
//                         port can have any value from [0~3]
//                         0      = Maximum Possible Instances that is 4 instance for above example
//                         1~3    = No of Instances to be calculated
//                         4~onwd = equivalent to Maximum no. of Instances
//  Output Port
//  1.Output Data Port   : Port is used for Data Output of given width, Order of Output Data is Dependent on value of 'ORDER' in template
//
//                                   ________________________
//                                  |                        |
//                                  |                        |
//    INPUT DATA------------------->|     DYNAMIC  FFT       |
//                                  |         WITH           |
//                                  |      GENERIC MIX       |
//    DYNAMIC REDUCTION------------>|         RADIX          |------------>OUTPUT DATA
//                                  |        IN-PLACE        |
//                                  |      ARCHITECTURE      |
//    INSTANCES COUNT-------------->|                        |
//                                  |                        |
//                                  |________________________|
//
//
//**************************************************************************************************************************

template < unsigned N_FFT, int RADIX, int ORDER, int TWIDDLE_PREC, int DIF_D_P, int DIF_D_I > 
class ac_fft_dif_r2pX_dyn_inpl
{
private:
  typedef ac_fixed < DIF_D_P, DIF_D_I, true > dif_fxp_data;
  typedef ac_complex < dif_fxp_data > dif_input, dif_output;
  typedef ac_int < ac::log2_ceil < ac::log2_ceil < N_FFT >::val - ac::log2_ceil < RADIX >::val + 1 >::val, 0 > dyn_port;
  typedef ac_int < ac::log2_ceil < N_FFT >::val - ac::log2_ceil < RADIX >::val, 0 > ins_port;

  ac_fft_dif_r2pX_dyn_inpl_core < N_FFT, RADIX, ORDER, TWIDDLE_PREC, DIF_D_P, DIF_D_I > fft;
  dif_input bank[RADIX][N_FFT / RADIX];

  //   Above line Instantiate RAM inside the design, RAM must mapped to a dual port Memory
  //   RAM will be spliced into n dual-port memories where n = 'RADIX'
  //   for Example
  //
  //   N_FFT (FFT Max) = 32, Radix = 4
  //
  //     Initial State of Bank
  //   Bank 1  Bank 2  Bank 3  Bank 4
  //     [X]     [X]     [X]     [X]
  //     [X]     [X]     [X]     [X]
  //     [X]     [X]     [X]     [X]
  //     [X]     [X]     [X]     [X]
  //     [X]     [X]     [X]     [X]
  //     [X]     [X]     [X]     [X]
  //     [X]     [X]     [X]     [X]
  //     [X]     [X]     [X]     [X]

  // coverAssert() used for basic template assert condition. This helps to validate if object
  // of this class created by the user has right set of parameters defined for it. Code will assert
  // during run-time if incorrect template values are used.

  void coverAssert () {
#ifdef ASSERT_ON
    AC_ASSERT((N_FFT == 2) || (N_FFT == 4) || (N_FFT == 8) || (N_FFT == 16) || (N_FFT == 32) || (N_FFT == 64) || (N_FFT == 128) || (N_FFT == 256) || (N_FFT == 512) || (N_FFT == 1024) || (N_FFT == 2048) || (N_FFT == 4096) || (N_FFT == 8192), "Number of FFT points not supported");
    AC_ASSERT((RADIX == 2) || (RADIX == 4) || (RADIX == 8) || (RADIX == 16) || (RADIX == 32) || (RADIX == 64) || (RADIX == 128), "Radix Value not supported.");
    AC_ASSERT(TWIDDLE_PREC <= 32, "Twiddle factor bitwidth greater than 32");
    AC_ASSERT(TWIDDLE_PREC >=  2, "Twiddle factor bitwidth lesser than 2");
    AC_ASSERT(DIF_D_P >= DIF_D_I, "Stage total bitwidth lesser than integer bitwidth.");
#endif
#ifdef COVER_ON
    cover (TWIDDLE_PREC <= 5);
#endif
  }

  template < int n >
  ac_int < n, 0 > shiftcr (ac_int < n, 0 > &input, ac_int < n, 0 > &shft) {
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

  template < int a >
  ac_int < a, 0 > bitrevint (ac_int < a, 0 > &Num) {
    ac_int < a, 0 > Num_br;

#pragma unroll yes
    BITREVERSE: for (ac_int < a + 1, 0 > itrator = 0; itrator < a; itrator++) {
      Num_br[a - 1 - itrator] = Num[itrator]; // bitreversal ac_int
    }
    return Num_br;
  }
public:

  ac_fft_dif_r2pX_dyn_inpl () { // Constructor
    ac::init_array < AC_VAL_DC > (&bank[0][0], N_FFT);
  }

#pragma hls_design interface
  void CCS_BLOCK(run) (ac_channel < dif_input > &inst, ac_channel < dif_output > &outst, ac_channel < dyn_port > &dyn_red_chn, ac_channel < ins_port > &ins_count) {
    coverAssert ();

    // Constant Computation depending on class template
    const int logN = ac::log2_ceil < N_FFT >::val;
    const int logRad = ac::log2_ceil < RADIX >::val;
    const int mixR = logN % logRad;
    const int Nstage = logN / logRad + 1;

    if (dyn_red_chn.available(1) & ins_count.available(1))
    {
      dyn_port dyn_red = dyn_red_chn.read ();      // Reading Dynamic Reduction Port
      ins_port instance_count = ins_count.read (); // Reading Instance Count Port

      if (inst.available (((instance_count == 0) ? (ac_int < logN + 1, false >) (N_FFT) : (ac_int < logN + 1, false >) ((N_FFT >> dyn_red) * instance_count))))
      {
        // Computing Dynamic characteristic of FFT
        ac_int < logN + 1, false > N_FFT_dyn = N_FFT >> dyn_red;
        ac_int < logN, 0 > logN_dyn = logN - dyn_red;
        ac_int < logN, 0 > Nstage_dyn = Nstage - (dyn_red >= mixR) - (dyn_red >= (mixR + logRad)) - (dyn_red >= (mixR + 2 * logRad)) - (dyn_red >= (mixR + 3 * logRad)) - (dyn_red >= (mixR + 4 * logRad)) - (dyn_red >= (mixR + 5 * logRad)) - (dyn_red >= (mixR + 6 * logRad)) - (dyn_red >= (mixR + 7 * logRad)) - (dyn_red >= (mixR + 8 * logRad)) - (dyn_red >= (mixR + 9 * logRad)) - (dyn_red >= (mixR + 10 * logRad)) - (dyn_red >= (mixR + 11 * logRad)) - (dyn_red >= (mixR + 12 * logRad));
        ac_int < logRad, 0 > RADIX_R_dyn = Nstage_dyn * logRad - logN_dyn;
        ac_int < logRad, 0 > mixR_dyn = (RADIX_R_dyn == 0) ? (ac_int < logRad, 0 >) 0 : (ac_int < logRad, 0 >) (logRad - RADIX_R_dyn);

        ac_int < logN - logRad + 1, 0 > fix_zero_wire_wdth = ((1 << (logN - dyn_red - logRad)) - 1);

        const ac_int < logN - logRad, 0 > masker_rad = (N_FFT - 1) >> logRad;

        ac_int < logN - logRad, 0 > j_mask = masker_rad >> dyn_red;
        ac_int < logN - logRad, 0 > block_counter = 0;

        // Loop Info:INPUT_DYN_ADDR_BLOCK_LOOP
        //
        // Input Stream will read in this loop and place in internal RAM
        //
        // 12 11 10 9 8 7 6 5 4 3 2 1 0--->>
        // >>-->>--Input-stream-->>-->>--->>
        // INPUT_DYN_ADDR_BLOCK_LOOP: for( Used to filling RAM  )
        //                         {
        //                             Sample place in RAM in Such a manner so that Mixing of Stages become conflict-free
        //                             (example for N_FFT=32 and RADIX=4 )
        //
        //                       for   dyn_red= 0  and ins_count= 1    (Dynamic Points FFT = 32, One instance )
        //                             B1       B2     B3      B4
        //                           a>[0]     [16]    [8]     [24]
        //                           b>[1]     [17]    [9]     [25]
        //       CASE-1              c>[2]     [18]    [10]    [26]
        //                           d>[3]     [19]    [11]    [27]
        //                           e>[4]     [20]    [12]    [28]
        //                           f>[5]     [21]    [13]    [29]
        //                             [6]     [22]    [14]    [30]
        //                             [7]     [23]    [15]    [31]
        //
        //                                        OR
        //
        //                       for   dyn_red= 1  and ins_count= 0   (Dynamic FFT Points =16 and Instance count = 2 )
        //                             B1       B2     B3       B4
        //                           a>[0]   e>[4]     [8]     [12]
        //                           b>[1]   f>[5]     [9]     [13]
        //                           c>[2]     [6]     [10]    [14]
        //       CASE-2              d>[3]     [7]     [11]    [15]
        //                             [0]     [4]     [8]     [12]
        //                             [1]     [5]     [9]     [13]
        //                             [2]     [6]     [10]    [14]
        //                             [3]     [7]     [11]    [15]
        //
        //                          }

#pragma hls_pipeline_init_interval 1
        INPUT_DYN_ADDR_BLOCK_LOOP: for (ac_int < logN + 1, false > blk = 0; blk < N_FFT; blk++) {
          ac_int < logRad, 0 > bank_temp, bank_add;
          bank_temp = blk.template slc < logRad > (logN - logRad - dyn_red);
          bank_add = shiftcr < logRad > (bank_temp, RADIX_R_dyn);
          bank[bank_add][(blk & j_mask) | ((blk >> logRad) & (~j_mask))] = inst.read ();

          if ((blk & ((1 << (logN - dyn_red)) - 1)) == ((N_FFT - 1) >> dyn_red)) {
            block_counter++;
          }

          if ((instance_count != 0) && (block_counter == instance_count)) { break; }
        }

        fft.fftDifR2pXDynInplCore (N_FFT_dyn, logN_dyn, bank, Nstage_dyn, dyn_red, RADIX_R_dyn, instance_count);

        ac_int < logN, 0 > out_add = 0;
        ac_int < logN - logRad + 1, 0 > mem_intr, mem_add;
        ac_int < logRad, 0 > m_no_msb = 0, bank_add;
        block_counter = 0;

        //  Loop Info: OUTPUT_DYN_NATURAL_ADDR_BLOCK_LOOP
        //
        //  Output of FFT will be kept in RAM blocks, Output stream sequence can be Natural or Bit-reversed
        //
        //
        //  OUTPUT_DYN_NATURAL_ADDR_BLOCK_LOOP: for (Used for Fetch Output Stream)
        //                          {
        //                              After FFT Computation FFT Result return to Mem Banks in arranged in specific Addressing,
        //                              to Fetch these Elements in Natural sequence using Reverse-address algo two Cases shown below
        //                              (example for N_FFT=32 and RADIX=4 )
        //
        //                        for   dyn_red= 0  and ins_count= 1    (Dynamic Points FFT = 32, One instance )
        //                              B1       B2     B3      B4
        //                              [0] >a  [21]    [10]    [31]
        //                              [4]     [17]    [14]    [27]
        //                              [8]     [29]    [2] >c  [23]           OUTPUT STREAM
        //          CASE-1              [12]    [25]    [6]     [19]          ---->>>>>   3 2 1 0
        //                              [16]    [5]     [26]    [15]
        //                              [20]    [1] >b  [30]    [11]
        //                              [24]    [13]    [18]    [7]
        //                              [28]    [9]     [22]    [3] >d
        //
        //
        //
        //                        for   dyn_red= 1  and ins_count= 0   (Dynamic FFT Points =16 and Instance count = 2 )
        //                              B1       B2     B3       B4
        //                              [0] >a  [5]     [10]    [15]
        //                              [4]     [1] >b  [14]    [11]
        //                              [8]     [13]    [2] >c  [7]            OUTPUT STREAM
        //                              [12]    [9]     [6]     [3] >d        ---->>>>>   3 2 1 0
        //           CASE-2             [0]     [5]     [10]    [15]
        //                              [4]     [1]     [14]    [11]
        //                              [8]     [13]    [2]     [7]
        //                              [12]    [9]     [6]     [3]
        //
        //                           }
        //

        if (ORDER == 1) {
#pragma hls_pipeline_init_interval 1
          OUTPUT_DYN_NATURAL_ADDR_BLOCK_LOOP: for (ac_int < logN + 1, false > blk = 0; blk < N_FFT; blk++) {
            ac_int < logN + 1, false > i = blk & j_mask;
            m_no_msb = blk.template slc < logRad > (logN - logRad - dyn_red);
            ac_int < logN - logRad + 1, 0 > blk_r = ((blk >> logRad) & (~j_mask));
            bank_add = bitrevint < logRad > (m_no_msb);
            m_no_msb = bank_add;
#pragma unroll yes
            OUTPUT_LOOP_NATURAL_SLICE: for (int p = 0; p <= (logN - logRad); p = (p + logRad)) {
              out_add.set_slc (p, m_no_msb);
            }
            mem_intr = i;
            mem_intr = bitrevint < logN - logRad + 1 > (mem_intr);
            mem_intr = mem_intr >> (dyn_red + 1);
            mem_add = ((fix_zero_wire_wdth) & (mem_intr ^ out_add));
            mem_add = mem_add | blk_r;
            outst.write (bank[bank_add][mem_add]);

            if ((blk & ((1 << (logN - dyn_red)) - 1)) == ((N_FFT - 1) >> dyn_red)) {
              block_counter++;
            }
            if (((instance_count != 0) && (block_counter == instance_count)) || (blk == N_FFT - 1)) { break; }
          }
        }

        //  Loop Info:OUTPUT_DYN_BITREVERSE_ADDR_BLOCK_LOOP
        //
        //  Output of FFT will be kept in RAM blocks, Output stream sequence can be Natural or Bit-reversed
        //
        //
        //  OUTPUT_DYN_BITREVERSE_ADDR_BLOCK_LOOP: for (Used for Fetch Output Stream)
        //                          {
        //                              After FFT Computation FFT Result return to Mem Banks in arranged in specific Addressing,
        //                              to Fetch these Elements in bitreversed sequence using Reverse-address algo two Cases shown below
        //                              (example for N_FFT=32 and RADIX=4 )
        //
        //                        for   dyn_red= 0  and ins_count= 1    (Dynamic Points FFT = 32, One instance )
        //                              B1       B2     B3      B4
        //                              [0] >a  [21]    [10]    [31]
        //                              [4]     [17]    [14]    [27]
        //                              [8] >c  [29]    [2]     [23]
        //          CASE-1              [12]    [25]    [6]     [19]          ---->>>>>   24 8 16 0
        //                              [16]>b  [5]     [26]    [15]
        //                              [20]    [1]     [30]    [11]
        //                              [24]>d  [13]    [18]    [7]
        //                              [28]    [9]     [22]    [3]
        //
        //
        //
        //                        for   dyn_red= 1  and ins_count= 0   (Dynamic FFT Points =16 and Instance count = 2 )
        //                              B1       B2     B3       B4
        //                              [0] >a  [5]     [10]    [15]
        //                              [4] >c  [1]     [14]    [11]
        //                              [8] >b  [13]    [2]     [7]
        //                              [12]>d  [9]     [6]     [3]           ---->>>>>   12 4 8 0
        //           CASE-2             [0]     [5]     [10]    [15]
        //                              [4]     [1]     [14]    [11]
        //                              [8]     [13]    [2]     [7]
        //                              [12]    [9]     [6]     [3]
        //
        //                           }

        else if (ORDER == 0) {
#pragma hls_pipeline_init_interval 1
          OUTPUT_DYN_BITREVERSE_ADDR_BLOCK_LOOP: for (ac_int < logN + 1, false > blk = 0; blk < N_FFT; blk++) {
            ac_int < logN + 1, false > i = ((blk >> logRad) & j_mask);
            m_no_msb = blk.template slc < logRad > (0);
            ac_int < logN - logRad + 1, 0 > blk_r = ((blk >> logRad) & (~j_mask));
#pragma unroll yes
            OUTPUT_LOOP_BITREVERSE_SLICE: for (int p = 0; p <= (logN - logRad); p = (p + logRad)) {
              out_add.set_slc (p, m_no_msb);
            }
            mem_add = ((fix_zero_wire_wdth) & (i ^ out_add));
            mem_add = mem_add | blk_r;
            bank_add = m_no_msb;
            outst.write (bank[bank_add][mem_add]);

            if ((blk & ((1 << (logN - dyn_red)) - 1)) == ((N_FFT - 1) >> dyn_red)) {
              block_counter++;
            }
            if (((instance_count != 0) && (block_counter == instance_count)) || (blk == N_FFT - 1)) {
              break;
            }
          }
        }
      }
    }

  }

};

#endif

