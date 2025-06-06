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
#ifndef __AC_CIC_DEC_FULL_PARAM_H
#define __AC_CIC_DEC_FULL_PARAM_H
#include <ac_fixed.h>
const unsigned R_TB = 7;
const unsigned M_TB = 2;
const unsigned N_TB = 4;

const int  in_W = 32;
const int  in_I = 16;
const bool in_S = true;

const int  out_W = 48;
const int  out_I = 32;
const bool out_S = true;

typedef ac_fixed <in_W, in_I, in_S> IN_TYPE_TB;
typedef ac_fixed <out_W, out_I, out_S> OUT_TYPE_TB;

#endif
