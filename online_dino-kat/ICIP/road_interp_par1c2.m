function [vars, opts] = road_interp_par1c2()
% Syntax:   [vars, opts] = road_interp_par1c2();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_interp_par1c2/data.mat';
vars.outpath  = 'road_interp_par1c2.mat';
vars.idim     = [128, nan, 200];
vars.SNR      = 19;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];

% interpVideo opts
opts.algo     = 'paint2';
opts.method   = 1;
