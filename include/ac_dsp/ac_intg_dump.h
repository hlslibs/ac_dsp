/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) DSP Library                                        *
 *                                                                        *
 *  Software Version: 3.8                                                 *
 *                                                                        *
 *  Release Date    : Tue May 13 15:37:10 PDT 2025                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.8.1                                               *
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
/***************************************************************************************
// File:        ac_intg_dump.h
//
// Description:
// The filter is designed by using class based approach. This file is C wrapper over
// core functionality of the filter. Member function of the class arun
// is the top function for the C++ design.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_intg_dump.h>
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Define input, output, sample type
//      typedef ac_fixed<32, 16, true> IN_TYPE;
//      typedef ac_fixed<64, 32, true> OUT_TYPE;
//      typedef ac_int < 16, false > N_TYPE;
//      const int NS  = 1024;
//      const int CHN = 4;
//      ac_intg_dump< IN_TYPE, ACC_TYPE, OUT_TYPE, N_TYPE, NS, CHN > filter_design1;
//      // Initialize channels for input, output, samples
//      ac_channel<IN_TYPE>  input;
//      ac_channel<OUT_TYPE> output;
//      ac_channel<N_TYPE> n_sample;
//      // Make sure to initialize the inputs and the samples
//      // Call the top-level function.
//      filter_design1.run(input, output, n_sample);
//
//      CCS_RETURN(0);
//    }
//
/*******************************************************************************************/

#ifndef _INCLUDED_AC_INTG_DUMP_H_
#define _INCLUDED_AC_INTG_DUMP_H_

#include <ac_channel.h>

#include <mc_scverify.h>
/*
 * class Declarations
 */

//===============================================================================================================
// Class: ac_intg_dump_core
//---------------------------------------------------------------------------------------------------------------
template < class IN_TYPE, class ACC_TYPE, class OUT_TYPE, class N_TYPE, int NS, int CHN >
class ac_intg_dump_core
{
private: // Data
  ACC_TYPE temp[CHN];

public: // Functions
  // Constructor
  ac_intg_dump_core() {
    for (int i = 0; i < CHN; i++) {
      temp[i] = 0.0; // Temporary workaround to avoid using ac::init_array and triggering a CDesignChecker error.
    }
  }

//------------------------------------------------------------------------------------------------------
// Member Function: intgDumpCore()

  void intgDumpCore(IN_TYPE data_in, OUT_TYPE &data_out, N_TYPE n_sample, int i, int j, bool &flag) {
    temp[i] = temp[i] + data_in;
    if (j == n_sample) {
      data_out = temp[i];
      temp[i] = 0.0;
      flag = true;
    }
  }
};
 
//===================================================================================================================
// Class: ac_intg_dump
// Description: This class contains the top functions which further class core functionality
//-------------------------------------------------------------------------------------------------------------------

template < class IN_TYPE, class ACC_TYPE, class OUT_TYPE, class N_TYPE, int NS, int CHN >
class ac_intg_dump
{
private: // Data
  ac_intg_dump_core < IN_TYPE, ACC_TYPE, OUT_TYPE, N_TYPE, NS, CHN > intg_dump_inst;

public: // Functions
  // Constructor
  ac_intg_dump() {}

//------------------------------------------------------------------------------------------------------
// Member Function: run()
// Description: run() is top function for C++ module.

#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void CCS_BLOCK(run)(ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, ac_channel < N_TYPE > &n_sample) {
    OUT_TYPE data_out_t;
    bool flag = false;
    N_TYPE n_sample_t;

#ifndef __SYNTHESIS__
    while (data_in.available(1))
#endif
    {
      n_sample_t = n_sample.read();
      ACC_LOOP: for (int j = 1; j <= NS; j++) {         // NS is maximum number of samples that should be accumulated data_in a configuration
        CHN_LOOP: for (int i = 0; i < CHN; i++) {       // Number of channels configured
          IN_TYPE data_in_t = data_in.read();
          intg_dump_inst.intgDumpCore(data_in_t, data_out_t, n_sample_t, i, j, flag);
          if (flag) {
            data_out.write(data_out_t);
          }
        }
        if (flag) { break; }
      }
      flag = false;
    }
  }
};

#endif

