/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) DSP Library                                        *
 *                                                                        *
 *  Software Version: 3.4                                                 *
 *                                                                        *
 *  Release Date    : Thu Nov 17 21:43:31 PST 2022                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.4.5                                               *
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
//************************************************************************************************************
//
// File: ac_cic_dec_full.h
//
// Description:
// 	The ac_cic_dec_full class serves as a C++ interface to the core integrator and differentiator chain
// 	classes located in "ac_cic_full_core.h". The member function "run()" contained in the class is the top
// 	function for the C++ design and calls the differentiator and integrator functions as hierarchical
// 	blocks. Differentiator and integrator functions create core class objects of the differentiator and
// 	integrator chain, respectively. The core classes are generic and intended to work well even with
// 	SystemC implementations in addition to the supplied C++ implementation.
//
// Usage:
// 	A sample testbench and its implementation look like this:
// 	#include <ac_cic_dec_full.h>
//
// 	#include <mc_scverify.h>
//
// 	CCS_MAIN(int arg, char **argc)
// 	{
// 		// Define input and output type.
// 		typedef ac_fixed<32, 16, true> IN_TYPE_TB;
//  	typedef ac_fixed<48, 32, true> OUT_TYPE_TB;
//  	// Initialize object of filter class with above input and output types, as well as:
//  	// Rate change factor = 7
//  	// Differential delay = 2
//		// Number of Stages = 4
// 		ac_cic_dec_full<IN_TYPE_TB, OUT_TYPE_TB, 7, 2, 4> filter_design1;
//  	// Initialize channels for input and output
//  	ac_channel<IN_TYPE_TB>  input;
//  	ac_channel<OUT_TYPE_TB> output;
//  	// Write 4 inputs to the input port
//  	input.write(.125);
//  	input.write(.375);
//  	input.write(1.525);
//  	input.write(.125);
//  	// Call the top-level function.
//  	filter_design1.run(input, output);
//
//		CCS_RETURN(0);
//    }
//
// Notes:
// 	Attempting to call the function with a type that is not implemented will result
// 	in a compile error.
// 	Currently, the block only accepts ac_fixed inputs and outputs.
// 	The file needs "ac_cic_full_core.h" in order to function.
//
//*********************************************************************************************************

#ifndef _INCLUDED_AC_CIC_DEC_FULL_H_
#define _INCLUDED_AC_CIC_DEC_FULL_H_

#include <ac_channel.h>
#include <ac_dsp/ac_cic_full_core.h>

#include <mc_scverify.h>

//=========================================================================
// Struct: power
// Description: templatized struct for computing power of a number 
//-------------------------------------------------------------------------

template <int base, int expon>
struct power {
  enum { value = base * power<base, expon - 1>::value };
};

//=========================================================================
// Struct: power
// Description: templatized struct for computing power of a number with a 
// zero component 
//-------------------------------------------------------------------------

template <int base>
struct power<base, 0> {
  enum { value = 1 };
};

//===============================================================================================================
// Struct: find_inter_type_cic_dec
// Description: provides parameterized bitwidths to ensure a lossless 
// intermediate type for the CIC filter. W, I and S are word width, 
// integer width and signedness of the input. 
//---------------------------------------------------------------------------------------------------------------
template <class IN_TYPE, unsigned R, unsigned M, unsigned N>
struct find_inter_type_cic_dec {
  enum {
    W    = IN_TYPE::width,
    I    = IN_TYPE::i_width,
    S    = IN_TYPE::sign,
    outF = W - I,
    // Two ways to calculate output bitwidth exist:
    // (i) outW = ceil(log2((R*M)^N)) + W + (!S) and (ii) outW = N*ceil(log2((R*M))) + W + (!S)
    // (i) Gives the exact bitwidth for lossless precision, and (ii) Gives an approximate bitwidth that is
    //     at most N - 1 bits larger than the exact amount.
    // If your design values of R, M and N are large, then the computation of ceil(log2((R*M)^N)) could throw
    // an overflow error at compile time. In that case it is recommended that equation (ii) be used, because
    // (R*M) is highly unlikely to cause an overflow when passed as a template parameter to the log2_ceil struct.
    // Or you can enter the value in manually. Make sure that you use only one equation and comment out the other
    // to avoid a compile time error.
    outW = ac::log2_ceil<power<R, N>::value *power<M, N>::value>::val + W + int(!S),   // ...(i)
    //outW = N*ac::log2_ceil<R *M>::val + W + int(!S),                                 // ...(ii)
    outI = outW - outF
  };
  typedef ac_fixed<outW, outI, true> INT_TYPE;
};

//=============================================================================================================
// Class: ac_cic_dec_full
// Description: This class contains the run() function, as well as integrator and differentiator functions for
// interpolation and decimation respectively. When all the stages of the filters have same precision,
// HW utilization can be achieved by partially unrolling low rate differentiator loop to allow adders
// to be reused without degrading the throughput.
//-------------------------------------------------------------------------------------------------------------

template < class IN_TYPE, class OUT_TYPE, unsigned R_, unsigned M_, unsigned N_ >
class ac_cic_dec_full
{
public: 
  //--------------------------------------------------------------------------
  // Constructor
  ac_cic_dec_full() : intg_inst(true) {}

//------------------------------------------------------------------------------------------------------
// Member Function: run()
// Description: run() is top function for C++ module. Based on filter type 
// configured it instantiates the integrator and comb sections as 
// hierarchical blocks.

#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void CCS_BLOCK(run)(ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out) {
    decIntg(data_in, inf);      // Integrator chain for decimation
    decDiff(inf, data_out);     // differentiator chain for decimation
  }

private:
  // Find a lossless intermediate type to be used for computations in the core design.
  typedef class find_inter_type_cic_dec <IN_TYPE, R_, M_, N_>::INT_TYPE INT_TYPE;

  // Instantiate object of integrator core class. Keep in mind that the output will be full-precision. 
  // The intermediate type will automatically use full precision bitwidths too.
  ac_cic_full_core_intg < IN_TYPE, INT_TYPE, R_, M_, N_ > intg_inst;
  // Instantiate object of differentiator core class in a similar manner
  ac_cic_full_core_diff < INT_TYPE, INT_TYPE, R_, M_, N_ > diff_inst;
  ac_channel < INT_TYPE > inf;  // interface between integrator and differentiator chain

//------------------------------------------------------------------------------------------------------
// Member Function: decIntg()
// Description: decIntg() is integrator section for decimation filter and for C++ module.
// It creates object of integrator core class "cic_f_core_intg" and calls its member function
// decIntgCore 

#pragma hls_pipeline_init_interval 1
#pragma hls_design
  void decIntg(ac_channel < IN_TYPE > &data_in, ac_channel < INT_TYPE > &data_out) {
    IN_TYPE data_in_t;
    INT_TYPE data_out_t;
#ifndef __SYNTHESIS__
    while (data_in.available(1))
#endif
    {
      data_in_t = data_in.read();
      intg_inst.decIntgCore(data_in_t, data_out_t);
      if (intg_inst.dvalid) {
        data_out.write(data_out_t);
      }
    }
  }

//------------------------------------------------------------------------------------------------------
// Member Function: decIntg()
// Description: decDiff is comb/differentiator section for decimation filter and for C++ module.
// It creates object of differentiator core calss "cic_full_core_diff" and calls its member function
// decDiffCore().
#pragma hls_pipeline_init_interval 1
#pragma hls_design
  void decDiff(ac_channel < INT_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out) {
    INT_TYPE data_in_t, data_out_t;
    OUT_TYPE data_out_final;
#ifndef __SYNTHESIS__
    while (data_in.available(1))
#endif
    {
      data_in_t = data_in.read();
      diff_inst.decDiffCore(data_in_t, data_out_t);
      // Convert the output of the core design, which is at full precision, to the output precision and write that to the output channel.
      data_out_final = data_out_t;
      data_out.write(data_out_final);
    }
  }
};

#endif

