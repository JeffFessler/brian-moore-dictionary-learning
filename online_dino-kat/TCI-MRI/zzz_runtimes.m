%% Batch runtimes 
% PINCAT
% 128 × 128 × 50

%
% Batch DINO-KAT, r = 1
%
% nIters = 75, stride = 2 x 2
%
% avg. time per outer iteration: 100.11s
%
%Iter[01] cost[3050.80] nrmse[0.261] deltaX[1.025e-01] deltaD[1.102e+00] deltaB[0.000e+00] sparsity[17.06%] time[140.78s] 
%Iter[02] cost[2259.72] nrmse[0.240] deltaX[3.795e-02] deltaD[4.001e-01] deltaB[9.977e-02] sparsity[14.90%] time[121.28s] 
%Iter[03] cost[1980.45] nrmse[0.223] deltaX[2.969e-02] deltaD[2.445e-01] deltaB[6.488e-02] sparsity[13.40%] time[114.92s] 
%Iter[04] cost[1772.79] nrmse[0.209] deltaX[2.768e-02] deltaD[1.495e-01] deltaB[5.692e-02] sparsity[11.90%] time[108.36s] 
%Iter[05] cost[1582.14] nrmse[0.195] deltaX[2.710e-02] deltaD[1.057e-01] deltaB[5.441e-02] sparsity[10.41%] time[97.32s] 
%Iter[06] cost[1397.62] nrmse[0.183] deltaX[2.649e-02] deltaD[1.605e-01] deltaB[5.289e-02] sparsity[8.95%] time[89.37s] 
%Iter[07] cost[1228.91] nrmse[0.174] deltaX[2.471e-02] deltaD[1.166e-01] deltaB[4.987e-02] sparsity[7.66%] time[80.78s]
%Iter[08] cost[1092.26] nrmse[0.167] deltaX[2.118e-02] deltaD[1.533e-01] deltaB[4.462e-02] sparsity[6.64%] time[72.59s]
%Iter[09] cost[994.77] nrmse[0.162] deltaX[1.674e-02] deltaD[1.472e-01] deltaB[3.801e-02] sparsity[5.93%] time[75.63s]
%
% [140.78, 121.28, 114.92, 108.36, 97.32, 89.37, 80.78, 72.59, 75.63]
%
addpath('../ISMRM/batch')
fcn = @batch_pincat_dinokat3;

% Run batch method
batch_run(fcn,3);

%% Online runtimes
% PINCAT
% 128 × 128 × 50

%
% Online DINO-KAT, r = 1
%
% nItersi = 50, nIters = 15, stride = 2 x 2, T = 5, dt = 1 (46 minibatches)
%
%Iter[01] ccost[101.28] nrmse[0.229] deltaX[0.000e+00] deltaD[0.000e+00] deltaB[1.134e-01] sparsity[8.49%] time[0.85s] 
%Iter[02] ccost[98.99] nrmse[0.229] deltaX[0.000e+00] deltaD[0.000e+00] deltaB[3.784e-02] sparsity[8.26%] time[0.83s] 
%Iter[03] ccost[98.19] nrmse[0.229] deltaX[0.000e+00] deltaD[0.000e+00] deltaB[2.344e-02] sparsity[8.14%] time[0.83s] 
%Iter[04] ccost[97.77] nrmse[0.229] deltaX[0.000e+00] deltaD[0.000e+00] deltaB[1.714e-02] sparsity[8.07%] time[0.85s] 
%Iter[05] ccost[97.50] nrmse[0.229] deltaX[0.000e+00] deltaD[0.000e+00] deltaB[1.416e-02] sparsity[8.02%] time[0.85s] 
%Iter[06] ccost[88.30] nrmse[0.217] deltaX[4.708e-02] deltaD[1.192e-01] deltaB[1.613e-02] sparsity[7.95%] time[1.44s] 
%Iter[07] ccost[85.87] nrmse[0.211] deltaX[1.736e-02] deltaD[5.826e-02] deltaB[3.525e-02] sparsity[7.47%] time[1.39s] 
%Iter[08] ccost[83.98] nrmse[0.206] deltaX[1.547e-02] deltaD[4.231e-02] deltaB[3.385e-02] sparsity[6.99%] time[1.41s] 
%Iter[09] ccost[82.27] nrmse[0.201] deltaX[1.474e-02] deltaD[3.872e-02] deltaB[3.308e-02] sparsity[6.54%] time[1.37s] 
%Iter[10] ccost[80.70] nrmse[0.197] deltaX[1.415e-02] deltaD[3.845e-02] deltaB[3.178e-02] sparsity[6.12%] time[1.36s] 
%Iter[11] ccost[79.25] nrmse[0.193] deltaX[1.342e-02] deltaD[3.860e-02] deltaB[3.048e-02] sparsity[5.74%] time[1.31s] 
%Iter[12] ccost[78.00] nrmse[0.189] deltaX[1.249e-02] deltaD[3.931e-02] deltaB[2.869e-02] sparsity[5.40%] time[1.28s] 
%Iter[13] ccost[76.88] nrmse[0.185] deltaX[1.164e-02] deltaD[3.529e-02] deltaB[2.725e-02] sparsity[5.11%] time[1.29s] 
%Iter[14] ccost[75.92] nrmse[0.182] deltaX[1.091e-02] deltaD[3.515e-02] deltaB[2.525e-02] sparsity[4.86%] time[1.26s] 
%Iter[15] ccost[75.09] nrmse[0.179] deltaX[1.009e-02] deltaD[3.544e-02] deltaB[2.370e-02] sparsity[4.64%] time[1.24s] 
%
% [0.85, 0.83, 0.83, 0.85, 0.85, 1.44, 1.39, 1.41, 1.37, 1.36, 1.31, 1.28, 1.29, 1.26, 1.24]
%
% avg. time per minibatch: 1.171s
%   ==> 53.85s total per outer iteration
%
fcn = @pincat_onlineDls_par3d;

%
% Online DINO-KAT, unitary
%
% nItersi = 50, nIters = 15, stride = 2 x 2, T = 5, dt = 1 (46 minibatches)
%
%Iter[01] ccost[45.83] nrmse[0.293] deltaX[3.101e-02] deltaD[1.123e+00] deltaB[1.470e-01] sparsity[20.71%] time[0.29s] 
%Iter[02] ccost[45.17] nrmse[0.291] deltaX[1.181e-02] deltaD[1.107e+00] deltaB[2.209e-02] sparsity[20.23%] time[0.26s] 
%Iter[03] ccost[44.82] nrmse[0.288] deltaX[9.560e-03] deltaD[1.100e+00] deltaB[1.742e-02] sparsity[19.89%] time[0.27s] 
%Iter[04] ccost[44.57] nrmse[0.286] deltaX[8.584e-03] deltaD[1.092e+00] deltaB[1.563e-02] sparsity[19.61%] time[0.29s] 
%Iter[05] ccost[44.36] nrmse[0.284] deltaX[8.040e-03] deltaD[1.095e+00] deltaB[1.487e-02] sparsity[19.35%] time[0.29s] 
%Iter[06] ccost[44.16] nrmse[0.282] deltaX[7.771e-03] deltaD[1.103e+00] deltaB[1.444e-02] sparsity[19.10%] time[0.28s] 
%Iter[07] ccost[43.98] nrmse[0.280] deltaX[7.605e-03] deltaD[1.101e+00] deltaB[1.419e-02] sparsity[18.86%] time[0.29s] 
%Iter[08] ccost[43.80] nrmse[0.278] deltaX[7.535e-03] deltaD[1.099e+00] deltaB[1.426e-02] sparsity[18.62%] time[0.30s] 
%Iter[09] ccost[43.62] nrmse[0.275] deltaX[7.502e-03] deltaD[1.100e+00] deltaB[1.423e-02] sparsity[18.37%] time[0.31s] 
%Iter[10] ccost[43.44] nrmse[0.273] deltaX[7.543e-03] deltaD[1.091e+00] deltaB[1.442e-02] sparsity[18.11%] time[0.27s] 
%Iter[11] ccost[43.25] nrmse[0.271] deltaX[7.699e-03] deltaD[1.103e+00] deltaB[1.481e-02] sparsity[17.84%] time[0.30s] 
%Iter[12] ccost[43.06] nrmse[0.269] deltaX[7.827e-03] deltaD[1.100e+00] deltaB[1.511e-02] sparsity[17.55%] time[0.28s] 
%Iter[13] ccost[42.86] nrmse[0.266] deltaX[7.946e-03] deltaD[1.094e+00] deltaB[1.525e-02] sparsity[17.25%] time[0.27s] 
%Iter[14] ccost[42.65] nrmse[0.264] deltaX[8.022e-03] deltaD[1.101e+00] deltaB[1.534e-02] sparsity[16.96%] time[0.28s] 
%Iter[15] ccost[42.43] nrmse[0.262] deltaX[8.225e-03] deltaD[1.096e+00] deltaB[1.590e-02] sparsity[16.64%] time[0.28s] 
%
% [0.29, 0.26, 0.27, 0.29, 0.29, 0.28, 0.29, 0.30, 0.31, 0.27, 0.30, 0.28, 0.27, 0.28, 0.28]
% avg. time per minibatch: 0.284s
%   ==> 13.06s total per outer iteration
%
fcn = @pincat_onlineUnitaryDls_par6b;

% Run online method
onlineDls_run(fcn,3);

%% Runtimes final
% PINCAT
% 128 × 128 × 50

% Computer
% 2016 MacBook Pro
% 2.7 GHz Intel Core i7
% 16 GB 2133 MHz LPDDR3

% MATLAB
% R2017b (9.3.0.713579)

% batch, r = 1, 2 x 2 spatial stride, 75 iterations
% avg. time per outer iteration: 100.11s

% online, r = 1, 2 x 2 spatial stride, minibatch size = 5, temporal stride = 1, 15 iterations
% avg. time per minibatch: 1.171s
% ==> avg. time per outer iteration: 53.85s
% ==> speed-up over batch: 1.9x

% online, unitary, 2 x 2 spatial stride, minibatch size = 5, temporal stride = 1, 15 iterations
% avg. time per minibatch: 0.284s
% ==> avg. time per outer iteration: 13.06s
% ==> speed-up over batch: 7.7x
