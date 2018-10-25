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