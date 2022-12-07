%% nIters0 = 0

onlineDls_run(@pincat_onlineDls_parTEST,1);

%{
***** Online Dictionary Least-Squares *****
Iter[01] cost[183.91] nrmse[0.283] deltaX[1.003e-01] deltaD[1.168e+00] deltaB[1.755e-01] sparsity[15.23%] time[2.05s] 
Iter[02] cost[138.78] nrmse[0.266] deltaX[3.771e-02] deltaD[3.520e-01] deltaB[9.838e-02] sparsity[13.57%] time[1.86s] 
Iter[03] cost[124.94] nrmse[0.253] deltaX[2.770e-02] deltaD[1.573e-01] deltaB[5.946e-02] sparsity[12.46%] time[1.86s] 
Iter[04] cost[115.67] nrmse[0.241] deltaX[2.443e-02] deltaD[1.171e-01] deltaB[5.077e-02] sparsity[11.50%] time[1.64s] 
Iter[05] cost[107.96] nrmse[0.230] deltaX[2.273e-02] deltaD[1.029e-01] deltaB[4.694e-02] sparsity[10.64%] time[1.59s] 
Iter[06] cost[101.05] nrmse[0.220] deltaX[2.167e-02] deltaD[9.982e-02] deltaB[4.492e-02] sparsity[9.81%] time[1.55s] 
Iter[07] cost[94.70] nrmse[0.211] deltaX[2.072e-02] deltaD[9.531e-02] deltaB[4.324e-02] sparsity[9.04%] time[1.56s] 
Iter[08] cost[88.73] nrmse[0.202] deltaX[1.997e-02] deltaD[9.253e-02] deltaB[4.235e-02] sparsity[8.29%] time[1.44s] 
Iter[09] cost[83.16] nrmse[0.194] deltaX[1.902e-02] deltaD[9.280e-02] deltaB[4.103e-02] sparsity[7.59%] time[1.42s] 
Iter[10] cost[77.95] nrmse[0.187] deltaX[1.827e-02] deltaD[9.078e-02] deltaB[3.952e-02] sparsity[6.94%] time[1.39s] 
Iter[11] cost[73.21] nrmse[0.180] deltaX[1.700e-02] deltaD[9.261e-02] deltaB[3.740e-02] sparsity[6.37%] time[1.35s] 
Iter[12] cost[69.14] nrmse[0.175] deltaX[1.537e-02] deltaD[9.161e-02] deltaB[3.483e-02] sparsity[5.88%] time[1.31s] 
Iter[13] cost[65.68] nrmse[0.170] deltaX[1.386e-02] deltaD[8.918e-02] deltaB[3.230e-02] sparsity[5.47%] time[1.28s] 
Iter[14] cost[62.87] nrmse[0.166] deltaX[1.210e-02] deltaD[8.553e-02] deltaB[2.925e-02] sparsity[5.14%] time[1.62s] 
Iter[15] cost[60.65] nrmse[0.162] deltaX[1.044e-02] deltaD[7.807e-02] deltaB[2.615e-02] sparsity[4.88%] time[1.38s] 
Iter[16] cost[58.88] nrmse[0.160] deltaX[9.053e-03] deltaD[7.058e-02] deltaB[2.325e-02] sparsity[4.68%] time[1.54s] 
Iter[17] cost[57.45] nrmse[0.157] deltaX[8.107e-03] deltaD[6.305e-02] deltaB[2.145e-02] sparsity[4.51%] time[1.47s] 
Iter[18] cost[56.23] nrmse[0.155] deltaX[7.264e-03] deltaD[5.866e-02] deltaB[2.000e-02] sparsity[4.36%] time[1.30s] 
Iter[19] cost[55.21] nrmse[0.153] deltaX[6.395e-03] deltaD[5.106e-02] deltaB[1.831e-02] sparsity[4.24%] time[1.36s] 
Iter[20] cost[54.36] nrmse[0.152] deltaX[5.772e-03] deltaD[4.731e-02] deltaB[1.690e-02] sparsity[4.14%] time[1.37s] 
***** Online Dictionary Least-Squares *****
Iter[01] cost[84.00] nrmse[0.157] deltaX[1.203e-01] deltaD[2.661e-01] deltaB[1.059e-01] sparsity[4.66%] time[1.70s] 
Iter[02] cost[73.40] nrmse[0.156] deltaX[1.680e-02] deltaD[1.633e-01] deltaB[4.355e-02] sparsity[4.27%] time[1.69s] 
Iter[03] cost[70.44] nrmse[0.155] deltaX[8.896e-03] deltaD[6.192e-02] deltaB[2.653e-02] sparsity[4.12%] time[1.50s] 
Iter[04] cost[68.98] nrmse[0.154] deltaX[6.530e-03] deltaD[3.935e-02] deltaB[2.017e-02] sparsity[4.02%] time[1.55s] 
Iter[05] cost[68.02] nrmse[0.153] deltaX[5.499e-03] deltaD[3.347e-02] deltaB[1.747e-02] sparsity[3.94%] time[1.39s] 
Iter[06] cost[67.25] nrmse[0.153] deltaX[4.923e-03] deltaD[3.429e-02] deltaB[1.585e-02] sparsity[3.87%] time[1.19s] 
Iter[07] cost[66.62] nrmse[0.152] deltaX[4.443e-03] deltaD[2.924e-02] deltaB[1.478e-02] sparsity[3.81%] time[1.50s] 
Iter[08] cost[66.07] nrmse[0.151] deltaX[4.018e-03] deltaD[1.978e-02] deltaB[1.375e-02] sparsity[3.77%] time[1.28s] 
Iter[09] cost[65.62] nrmse[0.151] deltaX[3.757e-03] deltaD[1.571e-02] deltaB[1.263e-02] sparsity[3.72%] time[1.16s] 
Iter[10] cost[65.24] nrmse[0.150] deltaX[3.406e-03] deltaD[1.465e-02] deltaB[1.187e-02] sparsity[3.69%] time[1.17s] 
%}

%% nIters0 = 5

onlineDls_run(@pincat_onlineDls_parTEST,1);

%{
***** Online Dictionary Least-Squares *****
Iter[01] cost[212.38] nrmse[0.292] deltaX[8.533e-02] deltaD[0.000e+00] deltaB[1.473e-01] sparsity[18.72%] time[0.71s] 
Iter[02] cost[180.21] nrmse[0.260] deltaX[5.444e-02] deltaD[0.000e+00] deltaB[9.032e-02] sparsity[15.57%] time[0.63s] 
Iter[03] cost[154.66] nrmse[0.231] deltaX[5.105e-02] deltaD[0.000e+00] deltaB[8.241e-02] sparsity[12.77%] time[0.65s] 
Iter[04] cost[130.50] nrmse[0.205] deltaX[4.991e-02] deltaD[0.000e+00] deltaB[7.971e-02] sparsity[10.09%] time[0.70s] 
Iter[05] cost[109.62] nrmse[0.185] deltaX[4.583e-02] deltaD[0.000e+00] deltaB[7.371e-02] sparsity[7.81%] time[0.64s] 
Iter[06] cost[98.35] nrmse[0.173] deltaX[5.262e-02] deltaD[1.224e+00] deltaB[1.137e-01] sparsity[5.92%] time[1.39s] 
Iter[07] cost[73.34] nrmse[0.162] deltaX[3.157e-02] deltaD[3.413e-01] deltaB[7.115e-02] sparsity[4.89%] time[1.42s] 
Iter[08] cost[64.68] nrmse[0.155] deltaX[2.089e-02] deltaD[2.042e-01] deltaB[4.432e-02] sparsity[4.31%] time[1.36s] 
Iter[09] cost[59.91] nrmse[0.151] deltaX[1.549e-02] deltaD[1.412e-01] deltaB[3.397e-02] sparsity[3.96%] time[1.27s] 
Iter[10] cost[56.80] nrmse[0.147] deltaX[1.215e-02] deltaD[1.118e-01] deltaB[2.794e-02] sparsity[3.72%] time[1.24s] 
Iter[11] cost[54.63] nrmse[0.144] deltaX[9.904e-03] deltaD[8.947e-02] deltaB[2.368e-02] sparsity[3.56%] time[1.25s] 
Iter[12] cost[53.05] nrmse[0.142] deltaX[8.212e-03] deltaD[7.648e-02] deltaB[2.043e-02] sparsity[3.43%] time[1.24s] 
Iter[13] cost[51.83] nrmse[0.141] deltaX[7.086e-03] deltaD[7.127e-02] deltaB[1.825e-02] sparsity[3.33%] time[1.21s] 
Iter[14] cost[50.83] nrmse[0.139] deltaX[6.149e-03] deltaD[6.139e-02] deltaB[1.658e-02] sparsity[3.25%] time[1.22s] 
Iter[15] cost[50.03] nrmse[0.138] deltaX[5.393e-03] deltaD[5.065e-02] deltaB[1.503e-02] sparsity[3.18%] time[1.20s] 
Iter[16] cost[49.36] nrmse[0.138] deltaX[4.825e-03] deltaD[4.442e-02] deltaB[1.365e-02] sparsity[3.13%] time[1.27s] 
Iter[17] cost[48.76] nrmse[0.137] deltaX[4.466e-03] deltaD[4.146e-02] deltaB[1.316e-02] sparsity[3.08%] time[1.20s] 
Iter[18] cost[48.23] nrmse[0.136] deltaX[4.040e-03] deltaD[3.858e-02] deltaB[1.234e-02] sparsity[3.03%] time[1.21s] 
Iter[19] cost[47.75] nrmse[0.136] deltaX[3.788e-03] deltaD[3.428e-02] deltaB[1.180e-02] sparsity[2.99%] time[1.21s] 
Iter[20] cost[47.32] nrmse[0.135] deltaX[3.599e-03] deltaD[3.162e-02] deltaB[1.128e-02] sparsity[2.95%] time[1.20s] 
***** Online Dictionary Least-Squares *****
Iter[01] cost[77.39] nrmse[0.150] deltaX[1.306e-01] deltaD[0.000e+00] deltaB[1.283e-01] sparsity[4.06%] time[1.04s] 
Iter[02] cost[69.45] nrmse[0.150] deltaX[1.596e-02] deltaD[0.000e+00] deltaB[4.479e-02] sparsity[3.51%] time[1.03s] 
Iter[03] cost[67.08] nrmse[0.150] deltaX[8.068e-03] deltaD[0.000e+00] deltaB[2.698e-02] sparsity[3.33%] time[1.03s] 
Iter[04] cost[65.92] nrmse[0.150] deltaX[5.365e-03] deltaD[0.000e+00] deltaB[1.962e-02] sparsity[3.23%] time[1.02s] 
Iter[05] cost[65.19] nrmse[0.150] deltaX[4.250e-03] deltaD[0.000e+00] deltaB[1.577e-02] sparsity[3.16%] time[1.05s] 
Iter[06] cost[63.10] nrmse[0.147] deltaX[9.529e-03] deltaD[7.242e-02] deltaB[1.426e-02] sparsity[3.11%] time[1.30s] 
Iter[07] cost[62.24] nrmse[0.145] deltaX[5.632e-03] deltaD[2.888e-02] deltaB[1.324e-02] sparsity[3.07%] time[1.25s] 
Iter[08] cost[61.67] nrmse[0.144] deltaX[4.237e-03] deltaD[1.979e-02] deltaB[1.134e-02] sparsity[3.04%] time[1.30s] 
Iter[09] cost[61.26] nrmse[0.143] deltaX[3.500e-03] deltaD[1.560e-02] deltaB[1.016e-02] sparsity[3.01%] time[1.27s] 
Iter[10] cost[60.95] nrmse[0.143] deltaX[3.015e-03] deltaD[1.296e-02] deltaB[9.052e-03] sparsity[2.99%] time[1.25s] 
%}

%% nIters0 = 5, Xmode = 1, Dmode = 0

onlineDls_run(@pincat_onlineDls_parTEST,1);

%{
***** Online Dictionary Least-Squares *****
Iter[01] cost[212.38] nrmse[0.292] deltaX[8.533e-02] deltaD[0.000e+00] deltaB[1.473e-01] sparsity[18.72%] time[0.75s] 
Iter[02] cost[180.21] nrmse[0.260] deltaX[5.444e-02] deltaD[0.000e+00] deltaB[9.032e-02] sparsity[15.57%] time[0.70s] 
Iter[03] cost[154.66] nrmse[0.231] deltaX[5.105e-02] deltaD[0.000e+00] deltaB[8.241e-02] sparsity[12.77%] time[0.64s] 
Iter[04] cost[130.50] nrmse[0.205] deltaX[4.991e-02] deltaD[0.000e+00] deltaB[7.971e-02] sparsity[10.09%] time[0.72s] 
Iter[05] cost[109.62] nrmse[0.185] deltaX[4.583e-02] deltaD[0.000e+00] deltaB[7.371e-02] sparsity[7.81%] time[0.63s] 
Iter[06] cost[98.35] nrmse[0.173] deltaX[5.262e-02] deltaD[1.224e+00] deltaB[1.137e-01] sparsity[5.92%] time[1.31s] 
Iter[07] cost[73.34] nrmse[0.162] deltaX[3.157e-02] deltaD[3.413e-01] deltaB[7.115e-02] sparsity[4.89%] time[1.30s] 
Iter[08] cost[64.68] nrmse[0.155] deltaX[2.089e-02] deltaD[2.042e-01] deltaB[4.432e-02] sparsity[4.31%] time[1.22s] 
Iter[09] cost[59.91] nrmse[0.151] deltaX[1.549e-02] deltaD[1.412e-01] deltaB[3.397e-02] sparsity[3.96%] time[1.19s] 
Iter[10] cost[56.80] nrmse[0.147] deltaX[1.215e-02] deltaD[1.118e-01] deltaB[2.794e-02] sparsity[3.72%] time[1.19s] 
Iter[11] cost[54.63] nrmse[0.144] deltaX[9.904e-03] deltaD[8.947e-02] deltaB[2.368e-02] sparsity[3.56%] time[1.17s] 
Iter[12] cost[53.05] nrmse[0.142] deltaX[8.212e-03] deltaD[7.648e-02] deltaB[2.043e-02] sparsity[3.43%] time[1.20s] 
Iter[13] cost[51.83] nrmse[0.141] deltaX[7.086e-03] deltaD[7.127e-02] deltaB[1.825e-02] sparsity[3.33%] time[1.17s] 
Iter[14] cost[50.83] nrmse[0.139] deltaX[6.149e-03] deltaD[6.139e-02] deltaB[1.658e-02] sparsity[3.25%] time[1.15s] 
Iter[15] cost[50.03] nrmse[0.138] deltaX[5.393e-03] deltaD[5.065e-02] deltaB[1.503e-02] sparsity[3.18%] time[1.17s] 
Iter[16] cost[49.36] nrmse[0.138] deltaX[4.825e-03] deltaD[4.442e-02] deltaB[1.365e-02] sparsity[3.13%] time[1.15s] 
Iter[17] cost[48.76] nrmse[0.137] deltaX[4.466e-03] deltaD[4.146e-02] deltaB[1.316e-02] sparsity[3.08%] time[1.24s] 
Iter[18] cost[48.23] nrmse[0.136] deltaX[4.040e-03] deltaD[3.858e-02] deltaB[1.234e-02] sparsity[3.03%] time[1.24s] 
Iter[19] cost[47.75] nrmse[0.136] deltaX[3.788e-03] deltaD[3.428e-02] deltaB[1.180e-02] sparsity[2.99%] time[1.23s] 
Iter[20] cost[47.32] nrmse[0.135] deltaX[3.599e-03] deltaD[3.162e-02] deltaB[1.128e-02] sparsity[2.95%] time[1.18s] 
***** Online Dictionary Least-Squares *****
Iter[01] cost[72.07] nrmse[0.148] deltaX[6.826e-02] deltaD[0.000e+00] deltaB[1.201e-01] sparsity[3.64%] time[0.99s] 
Iter[02] cost[67.72] nrmse[0.150] deltaX[1.078e-02] deltaD[0.000e+00] deltaB[3.450e-02] sparsity[3.36%] time[1.01s] 
Iter[03] cost[66.16] nrmse[0.150] deltaX[6.209e-03] deltaD[0.000e+00] deltaB[2.201e-02] sparsity[3.24%] time[0.99s] 
Iter[04] cost[65.29] nrmse[0.150] deltaX[4.591e-03] deltaD[0.000e+00] deltaB[1.719e-02] sparsity[3.17%] time[1.01s] 
Iter[05] cost[64.70] nrmse[0.150] deltaX[3.603e-03] deltaD[0.000e+00] deltaB[1.399e-02] sparsity[3.12%] time[0.99s] 
Iter[06] cost[62.73] nrmse[0.147] deltaX[9.277e-03] deltaD[7.165e-02] deltaB[1.365e-02] sparsity[3.06%] time[1.19s] 
Iter[07] cost[61.91] nrmse[0.145] deltaX[5.422e-03] deltaD[2.867e-02] deltaB[1.246e-02] sparsity[3.03%] time[1.20s] 
Iter[08] cost[61.39] nrmse[0.144] deltaX[4.015e-03] deltaD[1.944e-02] deltaB[1.066e-02] sparsity[3.00%] time[1.20s] 
Iter[09] cost[61.02] nrmse[0.143] deltaX[3.315e-03] deltaD[1.531e-02] deltaB[9.458e-03] sparsity[2.98%] time[1.21s] 
Iter[10] cost[60.71] nrmse[0.142] deltaX[2.850e-03] deltaD[1.348e-02] deltaB[9.051e-03] sparsity[2.96%] time[1.28s] 
%}

%% nIters0 = 5, Xmode = 1 ,Dmode = 1

onlineDls_run(@pincat_onlineDls_parTEST,1);

%{
***** Online Dictionary Least-Squares *****
Iter[01] cost[212.38] nrmse[0.292] deltaX[8.533e-02] deltaD[0.000e+00] deltaB[1.473e-01] sparsity[18.72%] time[0.71s] 
Iter[02] cost[180.21] nrmse[0.260] deltaX[5.444e-02] deltaD[0.000e+00] deltaB[9.032e-02] sparsity[15.57%] time[0.67s] 
Iter[03] cost[154.66] nrmse[0.231] deltaX[5.105e-02] deltaD[0.000e+00] deltaB[8.241e-02] sparsity[12.77%] time[0.66s] 
Iter[04] cost[130.50] nrmse[0.205] deltaX[4.991e-02] deltaD[0.000e+00] deltaB[7.971e-02] sparsity[10.09%] time[0.68s] 
Iter[05] cost[109.62] nrmse[0.185] deltaX[4.583e-02] deltaD[0.000e+00] deltaB[7.371e-02] sparsity[7.81%] time[0.76s] 
Iter[06] cost[98.35] nrmse[0.173] deltaX[5.262e-02] deltaD[1.224e+00] deltaB[1.137e-01] sparsity[5.92%] time[1.56s] 
Iter[07] cost[73.34] nrmse[0.162] deltaX[3.157e-02] deltaD[3.413e-01] deltaB[7.115e-02] sparsity[4.89%] time[1.43s] 
Iter[08] cost[64.68] nrmse[0.155] deltaX[2.089e-02] deltaD[2.042e-01] deltaB[4.432e-02] sparsity[4.31%] time[1.33s] 
Iter[09] cost[59.91] nrmse[0.151] deltaX[1.549e-02] deltaD[1.412e-01] deltaB[3.397e-02] sparsity[3.96%] time[1.31s] 
Iter[10] cost[56.80] nrmse[0.147] deltaX[1.215e-02] deltaD[1.118e-01] deltaB[2.794e-02] sparsity[3.72%] time[1.37s] 
Iter[11] cost[54.63] nrmse[0.144] deltaX[9.904e-03] deltaD[8.947e-02] deltaB[2.368e-02] sparsity[3.56%] time[1.26s] 
Iter[12] cost[53.05] nrmse[0.142] deltaX[8.212e-03] deltaD[7.648e-02] deltaB[2.043e-02] sparsity[3.43%] time[1.25s] 
Iter[13] cost[51.83] nrmse[0.141] deltaX[7.086e-03] deltaD[7.127e-02] deltaB[1.825e-02] sparsity[3.33%] time[1.25s] 
Iter[14] cost[50.83] nrmse[0.139] deltaX[6.149e-03] deltaD[6.139e-02] deltaB[1.658e-02] sparsity[3.25%] time[1.26s] 
Iter[15] cost[50.03] nrmse[0.138] deltaX[5.393e-03] deltaD[5.065e-02] deltaB[1.503e-02] sparsity[3.18%] time[1.24s] 
Iter[16] cost[49.36] nrmse[0.138] deltaX[4.825e-03] deltaD[4.442e-02] deltaB[1.365e-02] sparsity[3.13%] time[1.25s] 
Iter[17] cost[48.76] nrmse[0.137] deltaX[4.466e-03] deltaD[4.146e-02] deltaB[1.316e-02] sparsity[3.08%] time[1.25s] 
Iter[18] cost[48.23] nrmse[0.136] deltaX[4.040e-03] deltaD[3.858e-02] deltaB[1.234e-02] sparsity[3.03%] time[1.23s] 
Iter[19] cost[47.75] nrmse[0.136] deltaX[3.788e-03] deltaD[3.428e-02] deltaB[1.180e-02] sparsity[2.99%] time[1.27s] 
Iter[20] cost[47.32] nrmse[0.135] deltaX[3.599e-03] deltaD[3.162e-02] deltaB[1.128e-02] sparsity[2.95%] time[1.25s] 
***** Online Dictionary Least-Squares *****
Iter[01] cost[202.44] nrmse[0.133] deltaX[5.394e-02] deltaD[0.000e+00] deltaB[1.842e-01] sparsity[4.71%] time[0.73s] 
Iter[02] cost[197.92] nrmse[0.134] deltaX[1.928e-02] deltaD[0.000e+00] deltaB[3.656e-02] sparsity[4.23%] time[0.69s] 
Iter[03] cost[196.12] nrmse[0.135] deltaX[1.317e-02] deltaD[0.000e+00] deltaB[2.480e-02] sparsity[3.99%] time[0.62s] 
Iter[04] cost[195.12] nrmse[0.137] deltaX[9.905e-03] deltaD[0.000e+00] deltaB[1.938e-02] sparsity[3.84%] time[0.58s] 
Iter[05] cost[194.50] nrmse[0.138] deltaX[7.887e-03] deltaD[0.000e+00] deltaB[1.630e-02] sparsity[3.73%] time[0.60s] 
Iter[06] cost[95.57] nrmse[0.143] deltaX[3.554e-02] deltaD[1.304e+00] deltaB[1.002e-01] sparsity[3.07%] time[1.13s] 
Iter[07] cost[75.48] nrmse[0.141] deltaX[1.675e-02] deltaD[2.670e-01] deltaB[5.171e-02] sparsity[2.97%] time[1.12s] 
Iter[08] cost[69.84] nrmse[0.140] deltaX[1.090e-02] deltaD[1.092e-01] deltaB[3.028e-02] sparsity[2.84%] time[1.19s] 
Iter[09] cost[67.25] nrmse[0.138] deltaX[7.931e-03] deltaD[5.304e-02] deltaB[2.214e-02] sparsity[2.77%] time[1.22s] 
Iter[10] cost[65.80] nrmse[0.138] deltaX[6.116e-03] deltaD[3.442e-02] deltaB[1.729e-02] sparsity[2.73%] time[1.24s] 
%}
