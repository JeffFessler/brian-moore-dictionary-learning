function [vars, opts] = invivo_RPCA_R8_par1b()
% Syntax:   [vars, opts] = invivo_RPCA_R8_par1b();

% Knobs
vars.inpath   = 'invivo_R8.mat';
vars.rawpath  = 'invivo_RPCA_R8_par1b/data.mat';
vars.outpath  = 'invivo_RPCA_R8_par1b.mat';
vars.r        = nan;
vars.lambdaL  = 0.01  * logspace(-2,2,11);
vars.lambdaS  = 0.001 * logspace(-2,2,11);
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
