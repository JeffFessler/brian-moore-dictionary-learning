function [vars, opts] = road_interp_par1b()
% Syntax:   [vars, opts] = road_interp_par1b();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_interp_par1b/data.mat';
vars.outpath  = 'road_interp_par1b.mat';
vars.idim     = [256, nan, 200];
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];

% interpVideo opts
opts.algo     = 'paint2';
opts.method   = 1;
