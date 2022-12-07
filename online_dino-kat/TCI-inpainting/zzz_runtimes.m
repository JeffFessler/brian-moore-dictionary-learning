%% Batch runtimes 
% coastguard (noisy)
% 144 × 172 x 150

%
% Batch DINO-KAT, r = 1
%
% nIters = 20, stride = 4 x 4
%
% avg. time per outer iteration: 23.31s
%
%Iter[01] cost[4436.80] nrmse[0.128] deltaX[8.572e-02] deltaD[1.080e+00] deltaB[0.000e+00] sparsity[5.32%] time[28.20s] 
%Iter[02] cost[2982.66] nrmse[0.122] deltaX[2.633e-02] deltaD[4.465e-01] deltaB[7.388e-02] sparsity[3.74%] time[23.85s] 
%Iter[03] cost[2676.22] nrmse[0.119] deltaX[1.326e-02] deltaD[1.858e-01] deltaB[3.946e-02] sparsity[3.23%] time[22.37s] 
%Iter[04] cost[2538.83] nrmse[0.118] deltaX[9.028e-03] deltaD[1.245e-01] deltaB[2.762e-02] sparsity[2.97%] time[21.98s] 
%Iter[05] cost[2456.53] nrmse[0.116] deltaX[6.968e-03] deltaD[1.014e-01] deltaB[2.190e-02] sparsity[2.81%] time[21.83s] 
%Iter[06] cost[2400.66] nrmse[0.115] deltaX[5.664e-03] deltaD[8.557e-02] deltaB[1.838e-02] sparsity[2.71%] time[21.65s]
%
fcn = @dls_npar4b;

%
% Batch DINO-KAT, r = 5
%
% nIters = 20, stride = 4 x 4
%
% avg. time per outer iteration: 23.79s
%
%Iter[01] cost[3638.03] nrmse[0.131] deltaX[7.472e-02] deltaD[9.037e-01] deltaB[0.000e+00] sparsity[5.06%] time[28.56s] 
%Iter[02] cost[2828.04] nrmse[0.124] deltaX[2.547e-02] deltaD[4.463e-01] deltaB[5.763e-02] sparsity[3.59%] time[23.95s] 
%Iter[03] cost[2565.15] nrmse[0.121] deltaX[1.462e-02] deltaD[2.661e-01] deltaB[3.555e-02] sparsity[3.05%] time[24.03s] 
%Iter[04] cost[2435.41] nrmse[0.118] deltaX[9.936e-03] deltaD[1.849e-01] deltaB[2.608e-02] sparsity[2.79%] time[23.92s] 
%Iter[05] cost[2357.49] nrmse[0.117] deltaX[7.507e-03] deltaD[1.407e-01] deltaB[2.080e-02] sparsity[2.63%] time[21.25s] 
%Iter[06] cost[2304.94] nrmse[0.115] deltaX[6.007e-03] deltaD[1.179e-01] deltaB[1.748e-02] sparsity[2.52%] time[21.05s] 
%
fcn = @dls_npar4a;

% Run batch method
dls_run(fcn,3);

%% Online runtimes
% coastguard (noisy)
% 144 × 172 x 150

%
% Online DINO-KAT, r = 1
%
% nItersi = 50, nIters = 10, stride = 2 x 2, T = 5
%
% avg. time per minibatch: 0.43s
%   ==> 63s total per outer iteration w/ temporal stride of 1
%   ==> 12.9s total per outer iteration w/ temporal stride of 5
%
%Iter[01] ccost[44.14] nrmse[0.085] deltaX[4.300e-02] deltaD[2.702e-01] deltaB[2.785e-02] sparsity[1.51%] time[0.43s] 
%Iter[02] ccost[42.54] nrmse[0.084] deltaX[9.126e-03] deltaD[1.211e-01] deltaB[1.949e-02] sparsity[1.33%] time[0.43s] 
%Iter[03] ccost[42.20] nrmse[0.084] deltaX[5.059e-03] deltaD[1.268e-01] deltaB[1.052e-02] sparsity[1.29%] time[0.42s] 
%Iter[04] ccost[42.04] nrmse[0.084] deltaX[3.479e-03] deltaD[5.221e-02] deltaB[7.213e-03] sparsity[1.27%] time[0.43s] 
%Iter[05] ccost[41.95] nrmse[0.084] deltaX[2.639e-03] deltaD[4.026e-02] deltaB[5.771e-03] sparsity[1.26%] time[0.41s] 
%Iter[06] ccost[41.89] nrmse[0.084] deltaX[1.788e-03] deltaD[2.852e-02] deltaB[4.563e-03] sparsity[1.25%] time[0.42s] 
%Iter[07] ccost[41.85] nrmse[0.084] deltaX[1.459e-03] deltaD[2.207e-02] deltaB[3.996e-03] sparsity[1.24%] time[0.44s] 
%Iter[08] ccost[41.82] nrmse[0.084] deltaX[1.123e-03] deltaD[1.918e-02] deltaB[3.262e-03] sparsity[1.24%] time[0.42s] 
%Iter[09] ccost[41.80] nrmse[0.084] deltaX[8.345e-04] deltaD[1.484e-02] deltaB[2.991e-03] sparsity[1.24%] time[0.42s] 
%Iter[10] ccost[41.78] nrmse[0.084] deltaX[8.859e-04] deltaD[1.527e-02] deltaB[3.125e-03] sparsity[1.23%] time[0.42s] 
%
% nItersi = 50, nIters = 10, stride = 4 x 4, T = 5
%
% avg. time per minibatch: 0.09s
%   ==> 13.10s total per outer iteration w/ temporal stride of 1
%   ==>  2.70s total per outer iteration w/ temporal stride of 5
%
%Iter[01] ccost[17.64] nrmse[0.100] deltaX[8.594e-02] deltaD[3.463e-01] deltaB[9.034e-02] sparsity[3.54%] time[0.10s] 
%Iter[02] ccost[15.93] nrmse[0.096] deltaX[1.808e-02] deltaD[1.898e-01] deltaB[4.012e-02] sparsity[2.84%] time[0.09s] 
%Iter[03] ccost[15.33] nrmse[0.093] deltaX[1.103e-02] deltaD[1.149e-01] deltaB[2.572e-02] sparsity[2.57%] time[0.09s] 
%Iter[04] ccost[15.01] nrmse[0.092] deltaX[8.058e-03] deltaD[5.943e-02] deltaB[1.956e-02] sparsity[2.41%] time[0.09s] 
%Iter[05] ccost[14.82] nrmse[0.091] deltaX[6.175e-03] deltaD[5.299e-02] deltaB[1.513e-02] sparsity[2.32%] time[0.09s] 
%Iter[06] ccost[14.68] nrmse[0.090] deltaX[5.074e-03] deltaD[4.799e-02] deltaB[1.240e-02] sparsity[2.26%] time[0.09s] 
%Iter[07] ccost[14.59] nrmse[0.090] deltaX[4.262e-03] deltaD[4.454e-02] deltaB[1.109e-02] sparsity[2.21%] time[0.09s] 
%Iter[08] ccost[14.51] nrmse[0.089] deltaX[3.691e-03] deltaD[3.346e-02] deltaB[9.965e-03] sparsity[2.17%] time[0.08s] 
%Iter[09] ccost[14.44] nrmse[0.089] deltaX[3.351e-03] deltaD[3.148e-02] deltaB[9.757e-03] sparsity[2.13%] time[0.09s] 
%Iter[10] ccost[14.39] nrmse[0.089] deltaX[2.760e-03] deltaD[2.761e-02] deltaB[7.849e-03] sparsity[2.11%] time[0.09s]
%
fcn = @onlineDls_npar4d;

%
% Online DINO-KAT, r = 5
%
% nItersi = 50, nIters = 10, stride = 2 x 2, T = 5
%
% avg. time per minibatch: 0.355s
%   ==> 51.83s total per outer iteration w/ temporal stride of 1
%   ==> 10.65s total per outer iteration w/ temporal stride of 5
%
%Iter[01] ccost[43.52] nrmse[0.085] deltaX[4.618e-02] deltaD[1.824e-01] deltaB[2.582e-02] sparsity[1.42%] time[0.35s] 
%Iter[02] ccost[42.79] nrmse[0.085] deltaX[6.683e-03] deltaD[1.094e-01] deltaB[1.655e-02] sparsity[1.31%] time[0.36s] 
%Iter[03] ccost[42.62] nrmse[0.085] deltaX[3.167e-03] deltaD[4.799e-02] deltaB[9.888e-03] sparsity[1.28%] time[0.35s] 
%Iter[04] ccost[42.55] nrmse[0.085] deltaX[2.138e-03] deltaD[2.690e-02] deltaB[6.963e-03] sparsity[1.26%] time[0.36s] 
%Iter[05] ccost[42.50] nrmse[0.085] deltaX[1.582e-03] deltaD[1.912e-02] deltaB[5.613e-03] sparsity[1.25%] time[0.36s] 
%Iter[06] ccost[42.47] nrmse[0.085] deltaX[1.392e-03] deltaD[1.672e-02] deltaB[4.746e-03] sparsity[1.25%] time[0.35s] 
%Iter[07] ccost[42.45] nrmse[0.085] deltaX[1.037e-03] deltaD[1.187e-02] deltaB[4.058e-03] sparsity[1.24%] time[0.35s] 
%Iter[08] ccost[42.43] nrmse[0.085] deltaX[8.963e-04] deltaD[1.377e-02] deltaB[3.378e-03] sparsity[1.24%] time[0.36s] 
%Iter[09] ccost[42.42] nrmse[0.085] deltaX[8.805e-04] deltaD[9.450e-03] deltaB[3.538e-03] sparsity[1.24%] time[0.35s] 
%Iter[10] ccost[42.40] nrmse[0.085] deltaX[7.147e-04] deltaD[7.407e-03] deltaB[3.040e-03] sparsity[1.23%] time[0.36s]
%
% nItersi = 50, nIters = 10, stride = 4 x 4, T = 5
%
% avg. time per minibatch: 0.075s
%   ==> 10.95s total per outer iteration w/ temporal stride of 1
%   ==>  2.25s total per outer iteration w/ temporal stride of 5
%
%Iter[01] ccost[14.51] nrmse[0.093] deltaX[4.823e-02] deltaD[2.558e-01] deltaB[3.907e-02] sparsity[1.78%] time[0.08s] 
%Iter[02] ccost[13.75] nrmse[0.092] deltaX[1.123e-02] deltaD[1.537e-01] deltaB[1.807e-02] sparsity[1.64%] time[0.08s] 
%Iter[03] ccost[13.58] nrmse[0.091] deltaX[5.583e-03] deltaD[6.784e-02] deltaB[9.483e-03] sparsity[1.61%] time[0.07s] 
%Iter[04] ccost[13.50] nrmse[0.091] deltaX[3.601e-03] deltaD[5.024e-02] deltaB[7.170e-03] sparsity[1.59%] time[0.08s] 
%Iter[05] ccost[13.45] nrmse[0.091] deltaX[2.923e-03] deltaD[3.316e-02] deltaB[6.152e-03] sparsity[1.58%] time[0.07s] 
%Iter[06] ccost[13.42] nrmse[0.091] deltaX[2.271e-03] deltaD[2.362e-02] deltaB[5.374e-03] sparsity[1.57%] time[0.08s] 
%Iter[07] ccost[13.40] nrmse[0.091] deltaX[1.885e-03] deltaD[1.865e-02] deltaB[4.314e-03] sparsity[1.56%] time[0.08s] 
%Iter[08] ccost[13.38] nrmse[0.091] deltaX[1.600e-03] deltaD[1.709e-02] deltaB[4.292e-03] sparsity[1.56%] time[0.08s] 
%Iter[09] ccost[13.36] nrmse[0.091] deltaX[1.380e-03] deltaD[1.333e-02] deltaB[3.624e-03] sparsity[1.55%] time[0.07s] 
%Iter[10] ccost[13.35] nrmse[0.091] deltaX[1.390e-03] deltaD[1.255e-02] deltaB[3.628e-03] sparsity[1.55%] time[0.08s] 
%
fcn = @onlineDls_npar4a;

% Online DINO-KAT, unitary
%
% nItersi = 50, nIters = 10, stride = 2 x 2, T = 5
%
% avg. time per minibatch: 0.117s
%   ==> 17.08s total per outer iteration w/ temporal stride of 1
%   ==>  3.51s total per outer iteration w/ temporal stride of 5
%
%Iter[01] ccost[34.75] nrmse[0.084] deltaX[3.147e-02] deltaD[3.746e-01] deltaB[5.227e-02] sparsity[5.16%] time[0.11s] 
%Iter[02] ccost[32.55] nrmse[0.082] deltaX[1.167e-02] deltaD[3.078e-01] deltaB[2.478e-02] sparsity[3.83%] time[0.12s] 
%Iter[03] ccost[31.98] nrmse[0.080] deltaX[6.785e-03] deltaD[2.129e-01] deltaB[1.305e-02] sparsity[3.49%] time[0.12s] 
%Iter[04] ccost[31.74] nrmse[0.080] deltaX[4.751e-03] deltaD[1.520e-01] deltaB[9.008e-03] sparsity[3.35%] time[0.13s] 
%Iter[05] ccost[31.60] nrmse[0.079] deltaX[3.676e-03] deltaD[1.207e-01] deltaB[7.131e-03] sparsity[3.26%] time[0.12s] 
%Iter[06] ccost[31.51] nrmse[0.079] deltaX[2.977e-03] deltaD[1.004e-01] deltaB[5.855e-03] sparsity[3.21%] time[0.11s] 
%Iter[07] ccost[31.45] nrmse[0.078] deltaX[2.502e-03] deltaD[9.238e-02] deltaB[4.968e-03] sparsity[3.18%] time[0.12s] 
%Iter[08] ccost[31.40] nrmse[0.078] deltaX[2.191e-03] deltaD[7.692e-02] deltaB[4.396e-03] sparsity[3.16%] time[0.12s] 
%Iter[09] ccost[31.37] nrmse[0.078] deltaX[1.978e-03] deltaD[7.560e-02] deltaB[3.943e-03] sparsity[3.14%] time[0.12s] 
%Iter[10] ccost[31.34] nrmse[0.078] deltaX[1.783e-03] deltaD[6.748e-02] deltaB[3.468e-03] sparsity[3.12%] time[0.12s] 
%
% nItersi = 50, nIters = 10, stride = 4 x 4, T = 5
%
% avg. time per minibatch: 0.04s
%   ==> 5.84s total per outer iteration w/ temporal stride of 1
%   ==> 1.20s total per outer iteration w/ temporal stride of 5
%
%Iter[01] ccost[11.08] nrmse[0.099] deltaX[3.131e-02] deltaD[2.036e-01] deltaB[7.741e-02] sparsity[7.63%] time[0.04s] 
%Iter[02] ccost[10.76] nrmse[0.098] deltaX[1.196e-02] deltaD[1.160e-01] deltaB[2.267e-02] sparsity[6.54%] time[0.04s] 
%Iter[03] ccost[10.65] nrmse[0.097] deltaX[7.489e-03] deltaD[7.042e-02] deltaB[1.406e-02] sparsity[6.13%] time[0.04s] 
%Iter[04] ccost[10.59] nrmse[0.096] deltaX[5.610e-03] deltaD[5.282e-02] deltaB[1.063e-02] sparsity[5.90%] time[0.04s] 
%Iter[05] ccost[10.55] nrmse[0.096] deltaX[4.444e-03] deltaD[4.292e-02] deltaB[8.514e-03] sparsity[5.75%] time[0.04s] 
%Iter[06] ccost[10.53] nrmse[0.095] deltaX[3.744e-03] deltaD[3.669e-02] deltaB[7.081e-03] sparsity[5.65%] time[0.04s] 
%Iter[07] ccost[10.51] nrmse[0.095] deltaX[3.256e-03] deltaD[3.177e-02] deltaB[6.113e-03] sparsity[5.57%] time[0.05s] 
%Iter[08] ccost[10.50] nrmse[0.095] deltaX[2.873e-03] deltaD[3.027e-02] deltaB[5.462e-03] sparsity[5.50%] time[0.04s] 
%Iter[09] ccost[10.49] nrmse[0.094] deltaX[2.593e-03] deltaD[2.705e-02] deltaB[4.936e-03] sparsity[5.45%] time[0.04s] 
%Iter[10] ccost[10.48] nrmse[0.094] deltaX[2.281e-03] deltaD[2.422e-02] deltaB[4.437e-03] sparsity[5.41%] time[0.04s] 
%
fcn = @onlineDls_npar4g;

% Run online method
onlineDls_run(fcn,3);

%% Runtimes final
% dataset: coastguard (noisy)
% data size: 144 × 172 x 150

% Computer
% 2016 MacBook Pro
% 2.7 GHz Intel Core i7
% 16 GB 2133 MHz LPDDR3

% MATLAB
% R2017b (9.3.0.713579)

% batch, r = 1, 4 x 4 spatial stride
% avg. time per outer iteration: 23.31s

% online, r = 1, 2 x 2 spatial stride, minibatch size = 5, temporal stride = 5
% avg. time per minibatch: 0.43s
% ==> avg. time per outer iteration: 12.90s
% ==> speed-up over batch: 1.8x

% online, r = 1, 4 x 4 spatial stride, minibatch size = 5, temporal stride = 5
% avg. time per minibatch: 0.09s
% ==> avg. time per outer iteration: 2.70s
% ==> speed-up over batch: 8.6x

% online, unitary, 2 x 2 spatial stride, minibatch size = 5, temporal stride = 5
% avg. time per minibatch: 0.12s
% ==> avg. time per outer iteration: 3.60s
% ==> speed-up over batch: 6.5x

% online, unitary, 4 x 4 spatial stride, minibatch size = 5, temporal stride = 5
% avg. time per minibatch: 0.04s
% ==> avg. time per outer iteration: 1.20s
% ==> speed-up over batch: 19.4x
