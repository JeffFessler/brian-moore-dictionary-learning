function [vars, opts] = interp_par3b()
% Syntax:   [vars, opts] = interp_par3b();

% Knobs
vars.inpath   = 'data/coastguard2.mat';
vars.rawpath  = 'interp_par3b/data.mat';
vars.outpath  = 'interp_par3b.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.SNR      = inf;

% interpVideo opts
opts.algo     = 'tri3';
opts.method   = 'natural';
