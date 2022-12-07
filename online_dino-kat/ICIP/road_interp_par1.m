function [vars, opts] = road_interp_par1()
% Syntax:   [vars, opts] = road_interp_par1();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_interp_par1/data.mat';
vars.outpath  = 'road_interp_par1.mat';
vars.idim     = [128, nan, 200];
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];

% interpVideo opts
opts.algo     = 'paint2';
opts.method   = 1;
