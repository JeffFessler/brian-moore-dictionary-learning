function [vars, opts] = road_rpca_par2()
% Syntax:   [vars, opts] = road_rpca_par2();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_rpca_par2/data.mat';
vars.outpath  = 'road_rpca_par2.mat';
vars.idim     = [128, nan, 200];
vars.nIters   = 100;
vars.algo     = 'paint2';
vars.method   = 1;
vars.SNR      = inf;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.lambda   = 0.0025;
vars.mu       = 1 / sqrt(128 * 171);

% robustPCA opts
opts.A        = 1;
%opts.M;
opts.T        = 1;
opts.r        = nan;
%opts.nIters
%opts.L0;
%opts.S0;
%opts.Xtrue
opts.accel    = true;
opts.tau      = 0.99 / (1 + opts.accel);
opts.flag     = 1;
