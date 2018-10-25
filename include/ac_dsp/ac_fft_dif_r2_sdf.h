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
// File: ac_fft_dif_r2_sdf.h
//
// Description:
//    Nomenclature:
//                           ac_fft_dif_r2_sdf
//                         /    /    |   \     \                           
//                        /   FFT    |  Radix-2 \                          
//                  C-view           |            Single Delay Feedback
//                        Decimation in Frequency
//
//    Organization:
//     The ac_fft_dif_r2_sdf class serves as a C++ interface to the core class. The
//     ac_fft_dif_r2_sdf_core class is generic and intended to work well even with SystemC
//     implementations in addition to the supplied C++ implementation. The core class
//     instantiates delay stages as an objects of the 'ac_fft_dif_r2_sdf_stage' class.
//     The stages are implemented as a cascaded Single Feedback structure.
//
//    Order of Input/Output:
//     Input -- Natural
//     Output-- Bit reversed
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_fft_dif_r2_sdf.h>
//
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Initialize object of FFT class with 4 FFT points, Memory threshold = 32,
//      // Twiddle bitwidth = 19, I/O bitwidth = 18, I/O Integer Width = 2.
//      ac_fft_dif_r2_sdf < 4, 32, 19, 18, 2 > fft_design1;
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
//      // Write 4 zero inputs to flush out the pipeline
//      input.write(IO_type(0, 0));
//      input.write(IO_type(0, 0));
//      input.write(IO_type(0, 0));
//      input.write(IO_type(0, 0));
//      // Call the top-level function once more, to flush out the pipeline
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

#ifndef _INCLUDED_AC_FFT_DIF_R2_SDF_H_
#define _INCLUDED_AC_FFT_DIF_R2_SDF_H_

#include <ac_int.h>
#include <ac_fixed.h>
#include <ac_complex.h>

//******************************************************************
// class "ac_fft_dif_r2_sdf_stage"
// class has member functions--
//
// butterfly( )   -- Radix-2 Butterfly Structure
// rescale( )     -- Scale by 1/2
// multRescale( ) -- Multiply with twiddle factor and scale by 1/2
// readMem( )     -- Read from Circular buffer or register
// writeMem( )    -- Write to Circular buffer or register
// stageRun( )    -- Executes FFT stage functionality
//******************************************************************

// Templatized ac_fft_dif_r2_sdf_stage  Class

template< int STAGE, int MEM_TH, class com_p, class complex_round, class com_rnd_ext, class com_mult_type, class complext, class complext_fix>
class ac_fft_dif_r2_sdf_stage
{
  ac_int<STAGE + 1, false> iterator;
  ac_int<STAGE + 1, false> addr;              // used to access memory address

  com_rnd_ext memory_shift[1 << STAGE];       // Memory implementation (must be mapped to RAM)
  com_rnd_ext shift_reg[1 << STAGE];          // Register implementation (must be mapped to registers)

  // Declaration of Member Functions

private:

  // function butterfly Compute Radix-2 Butterfly

  void butterfly(com_rnd_ext &x, com_rnd_ext &y) {
    com_rnd_ext temp_1, temp_2 ;
    temp_1 = x;
    temp_2 = y;
    x = (temp_1 + temp_2);         // Butterfly Computation
    y = (temp_1 - temp_2);
  }

  // function rescale hard-coded to scale down the data by 1/2

  complex_round rescale(const com_rnd_ext in) {
    complex_round tx ;

    // The function is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    tx.r() = (in.r() >> 1);
    tx.i() = (in.i() >> 1);

    return tx ;
  }

  // function 'multRescale' computes multiplication of twiddle factor and scales by 1/2

  complex_round multRescale(const com_rnd_ext x, const complext y) {
    com_mult_type temp_1, temp_2, tx ;
    complex_round out ;

    temp_1 = x;
    temp_2 = y;

    tx = temp_1 * temp_2;                   /* Twiddle Multiplication */

    // The output is hard-coded to scale by 1/2. However, the user can pass it as-is by
    // eliminating the right-shift in the following two equations.
    out.r() = (tx.r() >> 1);
    out.i() = (tx.i() >> 1);

    return out ;
  }

  // function 'readMem' reads data from shift Register or Circular buffer

  void readMem(com_rnd_ext &out) {
    if ((1 << STAGE)<(MEM_TH)) {
      out = shift_reg[(1 << STAGE) - 1];
    } else {
      out = memory_shift[addr & ((1 << STAGE) - 1)];
    }
  }

  // function 'writeMem' read data from shift Register or Circular buffer

  void writeMem(com_rnd_ext &in) {
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

public:
  // function 'stageRun' compute of stage computation
  void stageRun(com_p &x, const complext w[1 << STAGE]) {
    com_rnd_ext xt;
    com_rnd_ext yt;
    complex_round tmp_out ;
    com_rnd_ext tmp_write_d ;
    complext tw, twd;

    xt = x;
    readMem(yt);
    if ( iterator[STAGE] ) {
      butterfly(yt, xt);
      tmp_out = rescale(yt);
      tmp_write_d = xt;
    } else {
      ac_int<STAGE + 2, false> n;
      n = (iterator) & ((1 << STAGE) - 1);
      complext J = complext(0, 1);
      ac_int<STAGE + 1, false> t;
      // Extraction of twiddle value from 1/8 Cycle of Complex exponential
      t = (1&((n << 2) >> STAGE))? (ac_int<STAGE + 1, false> )((((1 << STAGE) >> 2)) - (n & ((((1 << STAGE) >> 2))-1))):(ac_int<STAGE + 1, false>)(n&((((1 << STAGE) >> 2)) - 1));
      twd = w[t];
      tw.r() = ((1 & ((n << 2) >> STAGE))|(1&((n << 1) >> STAGE)))?(complext_fix)(-twd.r()):(complext_fix)(twd.r()) ;
      tw.i() = ((!(1 & ((n << 2) >> STAGE)))|(1&((n << 1) >> STAGE)))?(complext_fix)(twd.i()):(complext_fix)(-twd.i()) ;
      tw = ((1 & ((n << 2) >> STAGE))^(1 & ((n << 1) >> STAGE))) && (STAGE >= 2)?(complext)(J * tw.conj()):(complext)tw;    /* Calculation of Twiddle factor value */
      tmp_out = multRescale(yt, tw);
      tmp_write_d = xt;
    }

    writeMem(tmp_write_d);

    x = tmp_out;

    iterator = iterator + 1;
    addr++;
    return ;
  }

  // Constructor of fft class
  // Reset Action will be done here

  ac_fft_dif_r2_sdf_stage () {
    iterator = 0;
    addr = 0;
    ac::init_array<AC_VAL_DC>(&memory_shift[0], 1 << STAGE);
  }

};

//***********************************************************************************
// Description:  class "ac_fft_dif_r2_sdf_core " can instantiate upto 12 stages
//
//               class has the member function--
//               fftDifR2SdfCore() -- Implements core functionality the FFT. For one
//                                    sample, all stages execute per call of function
//                                    'fftDifR2SdfCore'
//***********************************************************************************

template < unsigned N_FFT, int MEM_TH, int TWID_PREC, int DIF_D0_P, int DIF_D0_I >
class ac_fft_dif_r2_sdf_core
{
private:
  // Type definition for Multipliers, Accumulator and stage variable for all stage
  // based on template args
  typedef ac_fixed<TWID_PREC, 2, true, AC_RND_INF>               fix_point_tw;
  typedef ac_complex<fix_point_tw>                               complext;

  typedef ac_fixed<DIF_D0_P, DIF_D0_I, true>                     dif_fix_point;
  typedef ac_complex<dif_fix_point>                              dif_complex;
  typedef ac_fixed<DIF_D0_P, DIF_D0_I, true, AC_RND, AC_SAT>     dif_fix_round;
  typedef ac_complex<dif_fix_round>                              dif_complex_round;
  typedef ac_fixed<DIF_D0_P + 1, DIF_D0_I + 1, true>             dif_fix_round_ext;
  typedef ac_complex<dif_fix_round_ext>                          dif_complex_round_ext;
  typedef ac_fixed<DIF_D0_P + TWID_PREC - 1, 1 + DIF_D0_I, true> dif_fix_mul;
  typedef ac_complex<dif_fix_mul>                                dif_comp_mul;

  // creating stage objects for FFT

  ac_fft_dif_r2_sdf_stage <11, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_11;
  ac_fft_dif_r2_sdf_stage <10, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_10;
  ac_fft_dif_r2_sdf_stage < 9, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_9;
  ac_fft_dif_r2_sdf_stage < 8, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_8;
  ac_fft_dif_r2_sdf_stage < 7, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_7;
  ac_fft_dif_r2_sdf_stage < 6, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_6;
  ac_fft_dif_r2_sdf_stage < 5, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_5;
  ac_fft_dif_r2_sdf_stage < 4, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_4;
  ac_fft_dif_r2_sdf_stage < 3, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_3;
  ac_fft_dif_r2_sdf_stage < 2, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_2;
  ac_fft_dif_r2_sdf_stage < 1, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_1;
  ac_fft_dif_r2_sdf_stage < 0, MEM_TH, dif_complex, dif_complex_round, dif_complex_round_ext, dif_comp_mul, complext, fix_point_tw> stage_0;

public:

  // Function will instantiate all stages of FFT according to No. of FFT-Points and pass twiddle vector.
  void fftDifR2SdfCore(const dif_complex &input, dif_complex &output) {
    bool stage_init = false;

#include<twiddles_20bits.h>

    dif_complex temp_11, temp_10, temp_9, temp_8, temp_7, temp_6, temp_5, temp_4, temp_3, temp_2, temp_1, temp_0;

    // creating Stage Objects for FFT

    if (N_FFT >= 4096) {
      temp_11 =   (dif_complex_round)input;
      stage_init = true;
      stage_11.stageRun(temp_11, twiddle_11);
    }

    if (N_FFT >= 2048) {
      if (stage_init) { temp_10 = (dif_complex_round)temp_11; }
      else            { temp_10 = (dif_complex_round)input;   }
      stage_init = true;
      stage_10.stageRun(temp_10, twiddle_10);
    }

    if (N_FFT >= 1024) {
      if (stage_init) { temp_9 = (dif_complex_round)temp_10; }
      else            { temp_9 = (dif_complex_round)input;   }
      stage_init = true;
      stage_9.stageRun(temp_9, twiddle_9);
    }

    if (N_FFT >= 512) {
      if (stage_init) { temp_8 = (dif_complex_round)temp_9; }
      else            { temp_8 = (dif_complex_round)input;  }
      stage_init = true;
      stage_8.stageRun(temp_8, twiddle_8);
    }

    if (N_FFT >= 256) {
      if (stage_init) { temp_7 = (dif_complex_round)temp_8; }
      else            { temp_7 = (dif_complex_round)input;  }
      stage_init = true;
      stage_7.stageRun(temp_7, twiddle_7);
    }

    if (N_FFT >= 128) {
      if (stage_init) { temp_6 = (dif_complex_round)temp_7; }
      else            { temp_6 = (dif_complex_round)input;  }
      stage_init = true;
      stage_6.stageRun(temp_6, twiddle_6);
    }

    if (N_FFT >= 64) {
      if (stage_init) { temp_5 = (dif_complex_round)temp_6; }
      else            { temp_5 = (dif_complex_round)input;  }
      stage_init = true;
      stage_5.stageRun(temp_5, twiddle_5);
    }

    if (N_FFT >= 32) {
      if (stage_init) { temp_4 = (dif_complex_round)temp_5; }
      else            { temp_4 = (dif_complex_round)input;  }
      stage_init = true;
      stage_4.stageRun(temp_4, twiddle_4);
    }

    if (N_FFT >= 16) {
      if (stage_init) { temp_3 = (dif_complex_round)temp_4; }
      else            { temp_3 = (dif_complex_round)input;  }
      stage_init = true;
      stage_3.stageRun(temp_3, twiddle_3);
    }

    if (N_FFT >= 8) {
      if (stage_init) { temp_2 = (dif_complex_round)temp_3; }
      else            { temp_2 = (dif_complex_round)input;  }
      stage_init = true;
      stage_2.stageRun(temp_2, twiddle_2);
    }

    if (N_FFT >= 4) {
      if (stage_init) { temp_1 = (dif_complex_round)temp_2; }
      else            { temp_1 = (dif_complex_round)input;  }
      stage_init = true;
      stage_1.stageRun(temp_1, twiddle_1);
    }

    if (N_FFT >= 2) {
      if (stage_init) { temp_0 = (dif_complex_round)temp_1; }
      else            { temp_0 = (dif_complex_round)input;  }
      stage_init = true;
      stage_0.stageRun(temp_0, twiddle_0);
    }

    output = (dif_complex_round)temp_0;

    stage_init = false;
  }
};

#include <ac_channel.h>
#ifdef ASSERT_ON
#include <ac_assert.h>
#endif
#ifdef COVER_ON
#include <ac_assert.h>
#endif

//************************************************************************************************
//
// class "ac_fft_dif_r2_sdf"
//
// class has member functions--
// void run( )                  -- Implements FFT C++-View
// void coverAssert( )          -- Contains basic Asserts and cover conditions
//
//************************************************************************************************

template < unsigned N_FFT, int MEM_TH, int TWID_PREC, int DIF_D0_P, int DIF_D0_I >
class ac_fft_dif_r2_sdf
{
public:
  bool write_out;

  typedef ac_fixed<DIF_D0_P, DIF_D0_I, true> dif_fx0;
  typedef ac_complex<dif_fx0> comp_dif;

  comp_dif b;

  ac_fft_dif_r2_sdf_core < N_FFT, MEM_TH, TWID_PREC, DIF_D0_P, DIF_D0_I >  fft ;

  ac_fft_dif_r2_sdf() { /*constructor*/
    write_out = false;
    b = 0;
  };

  // coverAssert() used for basic template assert condition. This helps to validate if object
  // of this class created in user code has right set of parameters defined for it. Code will assert
  // during compile time if incorrect template values are used.

  void coverAssert() {
#ifdef ASSERT_ON
    AC_ASSERT((N_FFT == 2) || (N_FFT == 4) || (N_FFT == 8) || (N_FFT == 16) || (N_FFT == 32) || (N_FFT == 64) || (N_FFT == 128) \
              || (N_FFT == 256) || (N_FFT == 512) || (N_FFT == 1024) || (N_FFT == 2048) || (N_FFT == 4096), "N_FFT is not a power of two");
    AC_ASSERT(TWID_PREC <= 32, "Twiddle bitwidth greater than 32");
#endif
#ifdef COVER_ON
    cover(TWID_PREC <= 5);
#endif
  }

  // Function 'run()' can be called in user code by instantiating the object of the class with specific template parameters.

#pragma hls_design interface
#pragma hls_pipeline_init_interval 1
  void run(ac_channel<comp_dif> &x1, ac_channel<comp_dif> &y1) {
    coverAssert();
    comp_dif a, y;

#pragma hls_pipeline_init_interval 1
    sample_loop:for (int i = 0; i < N_FFT; i++) {
      a = x1.read();
      fft.fftDifR2SdfCore(a, y);             // Calling Core of FFT Design
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

};

#endif

