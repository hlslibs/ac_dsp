/**************************************************************************
 *                                                                        *
 *  Algorithmic C (tm) DSP Library                                        *
 *                                                                        *
 *  Software Version: 3.4                                                 *
 *                                                                        *
 *  Release Date    : Mon Feb  6 09:12:03 PST 2023                        *
 *  Release Type    : Production Release                                  *
 *  Release Build   : 3.4.6                                               *
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
//***********************************************************************************************************//
// File:         ac_poly_dec.h
// Description:
//    This filter is designed using class based hierarchy. The member function "run()" contained in the
//    ac_poly_dec class is the top function for the C++ design.
//
// Usage:
//    A sample testbench and its implementation look like this:
//    #include <ac_poly_dec.h>
//    #include <mc_scverify.h>
//
//    CCS_MAIN(int arg, char **argc)
//    {
//      // Define input, output, coefficient and accumulator type
//      typedef ac_fixed<32, 16, true> IN_TYPE;
//      typedef ac_fixed<64, 32, true> OUT_TYPE;
//      typedef ac_fixed<32, 16, true> COEFF_TYPE;
//      typedef ac_fixed<64, 32, true> ACC_TYPE;
//      typedef struct {
//        COEFF_TYPE coeffs[NTAPS * DF];
//      }
//      STR_COEFF_TYPE;
//      ac_poly_dec <IN_TYPE, COEFF_TYPE, STR_COEFF_TYPE, ACC_TYPE, OUT_TYPE, NTAPS, DF> filter_design1;
//      ac_channel<IN_TYPE>  input;
//      ac_channel<OUT_TYPE> output;
//      ac_channel<STR_COEFF_TYPE> coeffs_ch;
//      // Make sure to initialize the inputs and the coeffs
//      // Call the top-level function.
//      filter_design1.run(input, output, coeffs_ch);
//
//      CCS_RETURN(0);
//    }
// Revision History:
//    3.4.0  - Added CDesignChecker waivers/fixes for ac_dsp IP blocks.
//             Changes made in general:
//               - FXD violations were fixed by assigning fixed value
//                 initializations to ac_fixed variables.
//
//***********************************************************************************************************//

#ifndef _INCLUDED_AC_POLY_DEC_H_
#define _INCLUDED_AC_POLY_DEC_H_

#include <ac_fixed.h>
#include <ac_channel.h>
// Include header for CCS_BLOCK
#include <mc_scverify.h>

//===================================================================================================================
// Class: ac_poly_dec
// Description: This class contains the top function for the polyphase decimator filter ('run()'). 
//-------------------------------------------------------------------------------------------------------------------

template < class IN_TYPE, class COEFF_TYPE, class STR_COEFF_TYPE, class ACC_TYPE, class OUT_TYPE, int NTAPS, int DF >
class ac_poly_dec
{
public: // Functions
  // Constructor
  ac_poly_dec() {
    ac::init_array < AC_VAL_0 > (taps, NTAPS * DF);
    ac::init_array < AC_VAL_0 > (acc1, DF);
    acc = 0.0;
  };

//------------------------------------------------------------------------------------------------------
// Member Function: run()
// Description: run() is top function for C++ module.
//
#pragma hls_pipeline_init_interval 1
#pragma hls_design interface
  void CCS_BLOCK(run)( ac_channel < IN_TYPE > &data_in, ac_channel < OUT_TYPE > &data_out, ac_channel < STR_COEFF_TYPE > &coeffs_st ) {
#ifndef __SYNTHESIS__
    while (coeffs_st.available(1))
#endif
    {
      coeffs_t = coeffs_st.read();
    }

#ifndef __SYNTHESIS__
    while (data_in.available(DF))
#endif
    {
      DF: for (int df = DF - 1; df >= 0; df--) {
#pragma hls_unroll
        SHIFT: for (int i = NTAPS * DF - 1; i >= 0; i--) {
          taps[i] = (i == 0) ? data_in.read() : taps[i - 1];
        }
#pragma hls_unroll
        MAC: for (int tp = 0; tp < NTAPS; tp++) {
          int tap_index = tp + NTAPS * df;
          acc1[df] = acc1[df] + taps[tp * DF] * coeffs_t.coeffs[tap_index];
        }
        acc =  acc + acc1[df];
        acc1[df] = 0.0;
      }
      OUT_TYPE acc_t = acc;
      data_out.write(acc_t);
      acc = 0.0;
    }
  }

private: // Data
  IN_TYPE taps[NTAPS * DF];
  ACC_TYPE acc, acc1[DF];
  STR_COEFF_TYPE coeffs_t;

};

#endif

