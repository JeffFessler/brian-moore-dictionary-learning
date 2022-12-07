function [vars, opts] = road_interp_par1e()
% Syntax:   [vars, opts] = road_interp_par1e();

% Knobs
vars.inpath   = 'data/road100.mat';
vars.rawpath  = 'road_interp_par1e/data.mat';
vars.outpath  = 'road_interp_par1e.mat';
vars.idim     = [128, nan, 100];
vars.SNR      = 25;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];

% interpVideo opts
opts.algo     = 'paint2';
opts.method   = 1;
