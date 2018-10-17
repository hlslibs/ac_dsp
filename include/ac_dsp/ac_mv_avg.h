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

/************************************************************************************
// File:       ac_mv_avg.h
//
// Description:
// filter is designed by using class based approach. This file is C wrapper over
// core functionality of the filter. Member function of the class run()
// is the top function for the C++ design Core classes of the filter is in
// "ac_mv_avg_core.h" and can be used by C++ and SystemC wrapper.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_mv_avg.h>
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Define input, output, coefficient and accumulator type
//      typedef ac_fixed<32, 16, true> IN_TYPE;
//      typedef ac_fixed<64, 32, true> OUT_TYPE;
//      typedef ac_fixed < 16, 2, true > ACC_TYPE;
//      typedef ac_fixed<32, 16, true> COEFF_TYPE;
//      typedef ac_int < 10, false > S_TYPE;
//      const int NS  = 1024;
//      const int CHN = 4;
//      ac_mv_avg_wrapper < MAXSZ, TAPS, AC_TYPE, IN_TYPE, OUT_TYPE, ACC_TYPE,
//      COEFF_TYPE, S_TYPE > filter_design1;
//      // Initialize channels for input, output, samples
//      ac_channel<IN_TYPE>  input;
//      ac_channel<OUT_TYPE> output;
//      ac_channel<S_TYPE> n_sample;
//      // Make sure to initialize the inputs and the number of samples
//      // Call the top-level function.
//      filter_design1.run(input, output, n_sample);
//
//      CCS_RETURN(0);
//    }
//
*************************************************************************************/

#ifndef _INCLUDED_AC_MV_AVG_H_
#define _INCLUDED_AC_MV_AVG_H_

#include <ac_fixed.h>
#include <ac_window.h>
#include <ac_channel.h>

/*
 * Class declarations:
 * class ac_mv_avg_core is the core class of the filter
 */

template < int MAX_SAMPLE, int TAPS, ac_window_mode WIN_TYPE, class IN_TYPE, class OUT_TYPE, class ACC_TYPE, class COEFF_TYPE >
class ac_mv_avg_core
{

private: // Data
  ac_window_1d_flag < IN_TYPE, TAPS, WIN_TYPE > w;
  ACC_TYPE acc_reg;
  const COEFF_TYPE *const coeffs;

public: // Functions
  // Constructor
  ac_mv_avg_core(const COEFF_TYPE *const ptr) : coeffs(ptr) {
    acc_reg = 0;
  }

  // mvAvgCore()  is core function of the filter and it multiplies with weights and accumulates the input data
  void mvAvgCore(IN_TYPE data_in, ac_int < 1, false > sol, ac_int < 1, false > eol, OUT_TYPE &data_out, bool &valid) {
    w.write(data_in, sol, eol);
    valid = w.valid();
    if (w.valid()) {
#pragma unroll
      ACC: for (int j = -TAPS / 2; j <= TAPS / 2; j++) {
        acc_reg = acc_reg + (ACC_TYPE) w[j] * coeffs[j + (TAPS / 2)]; /* Accumulator */
      }
      data_out = acc_reg;
      acc_reg = 0;
    }
  }

};

/*
 * class Declarations
 */

/*
 * class "ac_mv_avg"
 */

template < int MAX_SAMPLE, int TAPS, ac_window_mode WIN_TYPE, class IN_TYPE, class OUT_TYPE, class ACC_TYPE, class COEFF_TYPE, class S_TYPE >
class ac_mv_avg
{

public:
  // This pointer is set to public so that the user can extract the coeffs array from the base class and use it for their own
  // purposes in their testbench.
  const COEFF_TYPE *const cff_ptr;

  // constructor with pointer of const coeff array as an arg
  ac_mv_avg(const COEFF_TYPE *const c_ptr) : cff_ptr(c_ptr) { }

#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  // run() is top function for C++
  void run(ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, ac_channel < S_TYPE > &n_sample) {
    ac_int < 1, false > sol;
    ac_int < 1, false > eol;
    bool valid;
    OUT_TYPE data_out_t;
    IN_TYPE data_in_t;
    S_TYPE sample;
    S_TYPE n_sample_t;
    ac_mv_avg_core < MAX_SAMPLE, TAPS, WIN_TYPE, IN_TYPE, OUT_TYPE, ACC_TYPE, COEFF_TYPE > mv_avg_inst(cff_ptr);    // Object of the Core class of the filter
#ifndef __SYNTHESIS__
    while (n_sample.available(1))
#endif
    {
      n_sample_t = n_sample.read();
      if (WIN_TYPE == AC_WIN) { // Does not have boundary condition
        sample = n_sample_t;
      } else {
        sample = n_sample_t + TAPS / 2;
      }                                                                      // For AC_CLIP and AC_MIRROR, has boundary condition
    }
#ifndef __SYNTHESIS__
    while (data_in.available(1))
#endif
    {
      SAMPLE:
      for (int cnt_sample = 0; cnt_sample < MAX_SAMPLE + TAPS / 2; cnt_sample++) {
        // Loop iterates for maximum number of data samples set in parameter
        // and breaks when number of samples as set in the port n_sample are read
        if (cnt_sample < n_sample_t) {
          data_in_t = data_in.read();
        }
        eol = (cnt_sample == n_sample_t - 1) ? 1 : 0;                       // End of the element
        sol = (cnt_sample == 0) ? 1 : 0;                                    // Start of the element
        mv_avg_inst.mvAvgCore(data_in_t, sol, eol, data_out_t, valid);      // Call core functionality of the filter
        if (valid) {
          data_out.write(data_out_t);
        }
        if (cnt_sample == sample -1) { break; }  // Break condition written in the end to avoid unnecessary loop iteration and clock cycle
      }
    }
  }

};

#endif

