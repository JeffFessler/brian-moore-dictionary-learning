function [vars, opts] = invivo_RPCA_R8_par3()
% Syntax:   [vars, opts] = invivo_RPCA_R8_par3();

% Knobs
vars.inpath   = 'invivo_R8.mat';
vars.rawpath  = 'invivo_RPCA_R8_par3/data.mat';
vars.outpath  = 'invivo_RPCA_R8_par3.mat';
vars.r        = nan;
vars.lambdaL  = 0.1 * logspace(-0.5,0.5,5);
vars.lambdaS  = 0.001 * logspace(-0.5,0.5,5);
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
