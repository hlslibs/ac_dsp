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

//**********************************************************************************************************
//
// File: ac_cic_full_core.h
//
// Description:
//    Implements integrator and differentiator chains for decimation or interpolation filter configurations.
//
// Notes:
//    This file is included in both the ac_cic_dec_full.h and ac_cic_intr_full.h files. It contains the core
//    integrator and differentiation chain classes for the computation of the CIC filter output. The core
//    classes are compatible not only with the supplied C++ implementation in the aforementioned header
//    files, but are also compatible with suitable SystemC implementations, hence making the core generic.
//
//**********************************************************************************************************

#ifndef _INCLUDED_AC_CIC_FULL_CORE_H_
#define _INCLUDED_AC_CIC_FULL_CORE_H_

#include <ac_fixed.h>

#ifndef __SYNTHESIS__
#include <iostream>
using namespace std;
#endif

//************************************************************************************
//
// Class: ac_cic_full_core_intg
// Description: Integrator chain class for full precision implementation
//
// Public Interface:
//   decIntgCore()  -- Integrator chain for decimation configuration of the filter
//   intrIntgCore() -- Integrator chain for interpolation configuration of the filter
//
//************************************************************************************

template < class IN_TYPE, class OUT_TYPE, unsigned R_, unsigned M_, unsigned N_ >
class ac_cic_full_core_intg
{
private:
  bool                 valid;
  ac_int < 8, false >  rate_cnt1;
  ac_int < 8, false >  rate_cnt;
  OUT_TYPE             intg_reg[N_];

  // Function: intStage
  // Description: Implements integrators for number of stages configured
  void intStage(OUT_TYPE data_in, OUT_TYPE &data_out) {
#pragma unroll yes
    INTG_STG: for (int i = N_ - 1; i > 0; i--) {
      intg_reg[i] = intg_reg[i] + intg_reg[i - 1];
    }
    intg_reg[0] = data_in + intg_reg[0];
    data_out = intg_reg[N_ - 1];
  }

public:
  bool              dvalid;
  ac_int<8, false>  cnt;

  // Constructor
  ac_cic_full_core_intg(bool value) {
    valid = true;
    rate_cnt = 0;
    dvalid = value;
    cnt = 0;
    rate_cnt1 = R_ - 1;
    // initialize intg_reg to 0
    ac::init_array < AC_VAL_0 > (intg_reg, N_);
  }

  // Function: decIntgCore
  // Description:
  //   Integrator chain for decimation configuration of the filter
  //   and calls intStage() for integrator chain.
  //   valid and dvalid signals insure R_ -1 data are are discarded
  void decIntgCore(IN_TYPE data_in, OUT_TYPE &data_out) {
    OUT_TYPE data_out_t;
    OUT_TYPE data_in_t;

    data_in_t = (OUT_TYPE) data_in;

    if (rate_cnt == 0) {
      valid = true;
    } else {
      valid = false;
    }

    intStage(data_in_t, data_out_t);

    if (valid) {
      data_out = data_out_t;
    }

    dvalid = valid;

    rate_cnt++;
    if (rate_cnt > R_ - 1) {
      rate_cnt = 0;
    }

  }

  // Function: intrIntgCore
  // Description:
  //   Integrator chain for interpolation configuration of the filter that calls intStage() for integrator chain.
  //   This function also ensures that the correct number of zeroes are inserted between samples.
  void intrIntgCore(OUT_TYPE data_in, OUT_TYPE &data_out) {
    OUT_TYPE data_in_t;

    if (rate_cnt1 == R_ - 1) {
      data_in_t = (OUT_TYPE) data_in;
      rate_cnt1 = 0;
      dvalid = false;
    } else if (rate_cnt1 == R_ - 2) {
      data_in_t = 0;
      rate_cnt1++;
      dvalid = true;
    } else {
      data_in_t = 0;
      rate_cnt1++;
      dvalid = false;
    }
    intStage(data_in_t, data_out);
  }
};

//**************************************************************************
//
// Class: ac_cic_full_core_diff
// Description: differentiator chain class for full precision implementation
//
// Public Interface:
//   decDiffCore() - Implements differentiator chain for decimation
//   intrDiffCore() - Implements differentiator chain for interpolation
//
//**************************************************************************

template < class IN_TYPE, class OUT_TYPE, unsigned R_, unsigned M_, unsigned N_ >
class ac_cic_full_core_diff
{
public: // Functions

  // Constructor
  ac_cic_full_core_diff() {
    // initialize comb_dly_ln to 0
    ac::init_array < AC_VAL_0 > (&comb_dly_ln[0][0], N_ * M_);
  }

  // Function: decDiffCore
  // Description:
  //   Implements differentiator chain for decimation configuration of the filter
  //   by calling comb() function
  void decDiffCore(OUT_TYPE data_in, OUT_TYPE &data_out) {
    OUT_TYPE data_in_t, data_out_t;
    data_in_t = data_in;
    comb(data_in_t, data_out_t);
    data_out = data_out_t;
  }

  // Function: intrDiffCore
  // Description:
  //   Implements differentiator chain for interpolation configuration of the filter
  //   by calling comb() function
  void intrDiffCore(IN_TYPE data_in, OUT_TYPE &data_out) {
    OUT_TYPE data_in_t, data_out_t;
    data_in_t = (OUT_TYPE) data_in;
    comb(data_in_t, data_out_t);
    data_out = data_out_t;
  }

private: // Data
  OUT_TYPE comb_dly_ln[N_][M_];

private: // Functions

  // Function: comb
  // Description:
  //   Implements multiple differentiator stages based on number of stages configured
  void comb(OUT_TYPE data_in, OUT_TYPE &data_out) {
    OUT_TYPE data_in_comb_stg;
    OUT_TYPE data_out_comb_stg;
#pragma unroll yes
    COMB: for (int i = 0; i < N_; i++) {
      if (i == 0) {
        data_in_comb_stg = data_in;           // for first diff stage input is the data input
      } else {
        data_in_comb_stg = data_out_comb_stg; // for other diff stage input is the data output from previous stages
      }
      diffStage(data_in_comb_stg, i, data_out_comb_stg);
    }
    data_out = data_out_comb_stg;
  }

  // Function: diffStage
  // Description: Single differentiator stage with differential delay M_
  void diffStage(OUT_TYPE data_in, ac_int < 8, false > k, OUT_TYPE &data_out) {
    data_out = data_in - comb_dly_ln[k][M_ - 1];
#pragma unroll yes
    COMB_B: for (int i = 0; i < M_; i++) {
      if (i != 0) {
        comb_dly_ln[k][i] = comb_dly_ln[k][i - 1];
      }
    }
    comb_dly_ln[k][0] = data_in;
  }

};

#endif

