clc;
close all;
clear variables;

inW = 32;
inI = 16;
R = 7;
M = 2;
N = 4;
Ninput = 10000
% Create a CIC decimator object
cicDecim = dsp.CICDecimator(R, M, N);

Fs   = 500;
Fdec = Fs/R;

frel1 = 0.05;
F1 = frel1*Fdec;
frel2 = 1;
F2 = frel2*Fdec;

% Make sure that the number of inputs are a multiple of R, in order to ensure 
% that the matlab testbench succeeds.
% The closest multiple of R that exceeds the number of input data points in
% the testbench is computed below.
Ninput = ceil(Ninput/R) * R;
%32*R;
% Obtain vector from 0 to Ninput - 1
n = 0 : Ninput - 1;
% Obtain lower and higher tone samples
sin1 = sin(2*pi*F1*n/Fs);
sin2 = sin(2*pi*F2*n/Fs);
% Compute lossless bitwidths for the CIC filter.
% ac::log2_ceil<power<R, N>::value*power<M, N>::value>::val + W + int(!S)
outW = ceil(log2((R*M)^N)) + inW;
outI = outW - (inW - inI);

% Find the max value that can be stored for the given input fixed point
% configuration (assuming a signed input)
quantum_fp = 2 ^ (inI - inW);
max_fp = (2 ^ (inI - 1)) - quantum_fp;
% Superimpose one tone on the other
sincombine = sin1 + sin2;
% Normalize to max_fp
max_sincombine = max(abs(sincombine));
sincombine = (sincombine / max_sincombine) * max_fp;

% Convert the input to a fixed point representation and pass it to the
% matlab CIC filter.
sincombine_fp = fi(sincombine, 1, inW, inW - inI);
output_CIC = cicDecim(sincombine_fp');

% Uncomment this section if you'd like to generate plots.
%{
fvtool(cicDecim, 'NormalizedFrequency', 'on');
title('Frequency response of CIC output');

figure()
plot(sin1);
title('Input 1');

figure()
plot(sin2);
title('Input 2');

figure()
plot(double(sincombine_fp));
title('Combined input in fixed point');

figure()
plot(double(output_CIC));
title('CIC output');
%}

% Write out the input and reference output to the relevant files. The fixed
% point values are converted to double as it is faster for Matlab to
% process and write double values as compared to fixed point ones.
dlmwrite('ac_cic_dec_full_input.txt', (double(sincombine_fp))', 'delimiter', '', 'precision', 128);
dlmwrite('ac_cic_dec_full_ref.txt', double(output_CIC), 'delimiter', '', 'precision', 128);

% Write configuration parameters to header
str_param = ...
['#ifndef __AC_CIC_DEC_FULL_PARAM_H\n', ...
 '#define __AC_CIC_DEC_FULL_PARAM_H\n', ...
 '#include <ac_fixed.h>\n', ...
 'const unsigned R_TB = ', num2str(R), ';\n', ...
 'const unsigned M_TB = ', num2str(M), ';\n', ...
 'const unsigned N_TB = ', num2str(N), ';\n', ...
 '\n', ...
 'const int  in_W = ', num2str(inW), ';\n', ...
 'const int  in_I = ', num2str(inI), ';\n', ...
 'const bool in_S = true;\n', ...
 '\n', ...
 'const int  out_W = ', num2str(outW), ';\n', ...
 'const int  out_I = ', num2str(outI), ';\n', ...
 'const bool out_S = true;\n', ...
 '\n', ...
 'typedef ac_fixed <in_W, in_I, in_S> IN_TYPE_TB;\n', ...
 'typedef ac_fixed <out_W, out_I, out_S> OUT_TYPE_TB;\n', ...
 '\n', ...
 '#endif'];

fID = fopen('ac_cic_dec_full_param.h', 'w');
fprintf(fID, str_param);
fclose(fID);

