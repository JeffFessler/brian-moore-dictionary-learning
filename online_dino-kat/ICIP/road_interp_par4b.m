function [vars, opts] = road_interp_par4b()
% Syntax:   [vars, opts] = road_interp_par4b();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_interp_par4b/data.mat';
vars.outpath  = 'road_interp_par4b.mat';
vars.idim     = [256, nan, 200];
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];

% interpVideo opts
opts.algo     = 'tri3';
opts.method   = 'natural';
