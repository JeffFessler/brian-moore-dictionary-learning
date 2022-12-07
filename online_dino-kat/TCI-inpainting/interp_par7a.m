function [vars, opts] = interp_par7a()
% Syntax:   [vars, opts] = interp_par7a();

% Knobs
vars.inpath   = 'data/gflower.mat';
vars.rawpath  = 'interp_par7a/data.mat';
vars.outpath  = 'interp_par7a.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.SNR      = inf;

% interpVideo opts
opts.algo     = 'grid2';
opts.method   = 'cubic';
