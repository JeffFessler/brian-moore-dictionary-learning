function [vars, opts] = road_interp_par4c()
% Syntax:   [vars, opts] = road_interp_par4c();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_interp_par4c/data.mat';
vars.outpath  = 'road_interp_par4c.mat';
vars.idim     = [128, nan, 200];
vars.SNR      = 25;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];

% interpVideo opts
opts.algo     = 'tri3';
opts.method   = 'natural';
