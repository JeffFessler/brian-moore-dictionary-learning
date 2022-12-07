function [vars, opts] = interp_par5a()
% Syntax:   [vars, opts] = interp_par5a();

% Knobs
vars.inpath   = 'data/gbus.mat';
vars.rawpath  = 'interp_par5a/data.mat';
vars.outpath  = 'interp_par5a.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.SNR      = inf;

% interpVideo opts
opts.algo     = 'grid2';
opts.method   = 'cubic';
