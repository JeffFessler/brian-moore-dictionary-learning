function [vars, opts] = interp_par7b()
% Syntax:   [vars, opts] = interp_par7b();

% Knobs
vars.inpath   = 'data/gflower.mat';
vars.rawpath  = 'interp_par7b/data.mat';
vars.outpath  = 'interp_par7b.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.SNR      = inf;

% interpVideo opts
opts.algo     = 'tri3';
opts.method   = 'natural';
