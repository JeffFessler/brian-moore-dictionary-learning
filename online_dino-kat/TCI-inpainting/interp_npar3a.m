function [vars, opts] = interp_npar3a()
% Syntax:   [vars, opts] = interp_npar3a();

% Knobs
vars.inpath   = 'data/coastguard2.mat';
vars.rawpath  = 'interp_npar3a/data.mat';
vars.outpath  = 'interp_npar3a.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.SNR      = 19; %% NOISY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpVideo opts
opts.algo     = 'grid2';
opts.method   = 'cubic';
