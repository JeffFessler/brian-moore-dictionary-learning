function [vars, opts] = road_interp_par2b()
% Syntax:   [vars, opts] = road_interp_par2b();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_interp_par2b/data.mat';
vars.outpath  = 'road_interp_par2b.mat';
vars.idim     = [256, nan, 200];
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];

% interpVideo opts
opts.algo     = 'grid2';
opts.method   = 'cubic';
