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
//
// File: ac_cic_intr_full.h
//
// Description:
//    The ac_cic_intr_full class serves as a C++ interface to the core integrator and differentiator chain
//    classes located in "ac_cic_full_core.h". The member function "run()" contained in the class is the top
//    function for the C++ design and calls the differentiator and integrator functions as hierarchical
//    blocks. Differentiator and integrator functions create core class objects of the differentiator and
//    integrator chain, respectively. The core classes are generic and intended to work well even with
//    SystemC implementations in addition to the supplied C++ implementation.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_cic_intr_full.h>
//
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Define input and output type.
//      typedef ac_fixed<32, 16, true> IN_TYPE_TB;
//      typedef ac_fixed<49, 33, true> OUT_TYPE_TB;
//      // Initialize object of filter class with above input and output types, as well as:
//      // Rate change factor = 7
//      // Differential delay = 2
//      // Number of Stages = 5
//      ac_cic_intr_full<IN_TYPE_TB, OUT_TYPE_TB, 7, 2, 5> filter_design1;
//      // Initialize channels for input and output
//      ac_channel<IN_TYPE_TB>  input;
//      ac_channel<OUT_TYPE_TB> output;
//      // Write 4 inputs to the input port
//      input.write(.125);
//      input.write(.375);
//      input.write(1.525);
//      input.write(.125);
//      // Call the top-level function.
//      filter_design1.run(input, output);
//
//      CCS_RETURN(0);
//    }
//
// Notes:
//    Attempting to call the function with a type that is not implemented will result
//    in a compile error.
//    Currently, the block only accepts ac_fixed inputs and outputs.
//    The file needs "ac_cic_full_core.h" in order to function.
//
//*********************************************************************************************************

#ifndef _INCLUDED_AC_CIC_INTR_FULL_H_
#define _INCLUDED_AC_CIC_INTR_FULL_H_

#include <ac_channel.h>
#include <ac_dsp/ac_cic_full_core.h>      /* Core functionality for the filter, common for C++ and System C */

// Templatized structs for computing power of a number.
template <int base, int expon>
struct power {
  enum { value = base * power<base, expon - 1>::value };
};

template <int base>
struct power<base, 0> {
  enum { value = 1 };
};

// This struct provides parameterized bitwidths to ensure a lossless intermediate type for the CIC filter.
// W, I and S are word width, integer width and signedness of the input.
template <class IN_TYPE, unsigned R, unsigned M, unsigned N>
struct find_inter_type_cic_intr {
  enum {
    W    = IN_TYPE::width,
    I    = IN_TYPE::i_width,
    S    = IN_TYPE::sign,
    outF = W - I,
    // Two ways to calculate output bitwidth exist:
    // (i) outW = ceil(log2((R^(N - 1)*(M^N)) + W + (!S) and (ii) outW = N*ceil(log2((R*M))) - floor(log2(R)) + W + (!S)
    // (i) Gives the exact bitwidth for lossless precision, and (ii) Gives an approximate bitwidth.
    // If your design values of R, M and N are large, then the computation of (R^(N - 1))*(M^N) could throw
    // an overflow error at compile time. In that case it is recommended that equation (ii) be used, because
    // (R*M) is highly unlikely to cause an overflow when passed as a template parameter to the log2_ceil struct.
    // Or you can enter the value in manually. If using an equation, make sure that you use only one equation and
    // comment out the other to avoid a compile time error.
    outW = ac::log2_ceil<power<R, N - 1>::value*power<M, N>::value>::val + W + int(!S),  // ...(i)
    //outW = N*ac::log2_ceil<R*M>::val - ac::log2_floor<R>::val + W + int(!S),           // ...(ii)
    outI = outW - outF
  };
  typedef ac_fixed<outW, outI, true> INT_TYPE;
};

//*******************************************************************************************************
//
//  class "ac_cic_intr_full"
//
//  This class contains the run() function, as well as, integrator and differentiator functions for
//  interpolation and decimation respectively. When all the stages of the filters have same precision,
//  HW utilization can be achieved by partially unrolling low rate differentiator loop to allow adders
//  to be reused without degrading the throughput.
//
//*******************************************************************************************************

template < class IN_TYPE, class OUT_TYPE, unsigned R_, unsigned M_, unsigned N_ >
class ac_cic_intr_full        /* CIC class for full precision implementation */
{
private:
  // Find a lossless intermediate type to be used for computations in the core design.
  typedef typename find_inter_type_cic_intr <IN_TYPE, R_, M_, N_>::INT_TYPE INT_TYPE;
  // Instantiate object of integrator core class. Keep in mind that the output will be full-precision. The intermediate type will automatically
  // use full precision bitwidths too.
  ac_cic_full_core_intg < IN_TYPE, INT_TYPE, R_, M_, N_ > intg_inst;
  // Instantiate object of differentiator core class in a similar manner
  ac_cic_full_core_diff < INT_TYPE, INT_TYPE, R_, M_, N_ > diff_inst;
  ac_channel < INT_TYPE > inf;   // interface between integrator and differentiator chain
#if !defined(__SYNTHESIS__) && defined(AC_CIC_INTR_FULL_H_DEBUG)
  bool print_once;
#endif

  // intrDiff is comb/differentiator section for interpolation filter and for C++ module.
  // It creates object of differentiator core class "cic_full_core_diff" and calls its member function
  // intrDiffCore. see ac_cic_full_core.h for intrDiffCore.
#pragma hls_pipeline_init_interval 1
#pragma hls_design
  void intrDiff(ac_channel < IN_TYPE > &data_in, ac_channel < INT_TYPE > &data_out) {
    IN_TYPE data_in_t;
    INT_TYPE data_out_t;

#ifndef __SYNTHESIS__
    while (data_in.available(1))
#endif
    {
      data_in_t = data_in.read();
      diff_inst.intrDiffCore(data_in_t, data_out_t);
      data_out.write(data_out_t);
    }
  }

  // intrIntg is integrator section for interpolation filter and for C++ module.
  // It creates object of inegrator core calss "cic_f_core_intg" and calls its member function
  // intrIntgCore(). see ac_cic_full_core.h for intrIntgCore().
#pragma hls_pipeline_init_interval 1
#pragma hls_design
  void intrIntg(ac_channel < INT_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out) {
    INT_TYPE data_in_t, data_out_t;
    OUT_TYPE data_out_final;

#ifndef __SYNTHESIS__
    while (data_in.available(1))
#endif
    {
      if (intg_inst.dvalid) {
        data_in_t = data_in.read();
      }
      intg_inst.intrIntgCore(data_in_t, data_out_t);
      // Convert the output of the core design, which is at full precision, to the output precision and write that to the output channel.
      data_out_final = data_out_t;
      if (intg_inst.cnt < N_ - 1) {
        intg_inst.cnt++;
      } else {
        data_out.write(data_out_final);
      }
    }
  }

public: // Functions
  // Constructor
  ac_cic_intr_full() : intg_inst(true) 
  {
#if !defined(__SYNTHESIS__) && defined(AC_CIC_INTR_FULL_H_DEBUG)
    print_once = true;
#endif
  }

  // run() is top function for C++ module. Based on filter type configured
  // it instantiates integrator and comb sections as hierarchical blocks.
#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void run(ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out) {
    intrDiff(data_in, inf);      // differentiator chain for Interpolation
    intrIntg(inf, data_out);     // Integrator chain for Interpolation

#if !defined(__SYNTHESIS__) && defined(AC_CIC_INTR_FULL_H_DEBUG)
    if (print_once) {
      print_once = false;
      cout << "INT_TYPE = " << INT_TYPE::type_name() << endl;
    }
#endif
  }

};

#endif

