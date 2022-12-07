function [vars, opts] = road_interp_par4e2()
% Syntax:   [vars, opts] = road_interp_par4e2();

% Knobs
vars.inpath   = 'data/road100.mat';
vars.rawpath  = 'road_interp_par4e2/data.mat';
vars.outpath  = 'road_interp_par4e2.mat';
vars.idim     = [128, nan, 100];
vars.SNR      = 19;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];

% interpVideo opts
opts.algo     = 'tri3';
opts.method   = 'natural';
