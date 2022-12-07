function [vars, opts] = pincat_RPCA_R8_par1n()
% Syntax:   [vars, opts] = pincat_RPCA_R8_par1n();

% Knobs
vars.inpath   = 'pincat_R8n.mat';
vars.rawpath  = 'pincat_RPCA_R8_par1n/data.mat';
vars.outpath  = 'pincat_RPCA_R8_par1n.mat';
vars.r        = nan;
vars.lambdaL  = 0.10 * logspace(-2,2,20);
vars.lambdaS  = 0.01 * logspace(-2,2,20);
vars.nIters   = 250;

% RPCA opts
%opts.A      = A;
%opts.T      = T;
%opts.r      = r;
%opts.nIters = vars.nIters;
%opts.L0     = reshape(data.Xfft,[],nt);
%opts.S0     = zeros(ny * nx,nt);
%opts.Xtrue  = reshape(data.Xtrue,[],nt);
opts.accel   = false;
opts.tau     = 1;
opts.flag    = 1;
