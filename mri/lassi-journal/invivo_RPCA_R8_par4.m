function [vars, opts] = invivo_RPCA_R8_par4()
% Syntax:   [vars, opts] = invivo_RPCA_R8_par4();

% Knobs
vars.inpath   = 'invivo_R8.mat';
vars.rawpath  = 'invivo_RPCA_R8_par4/data.mat';
vars.outpath  = 'invivo_RPCA_R8_par4.mat';
vars.r        = 1:10;
vars.lambdaL  = nan;
vars.lambdaS  = 0.001 * logspace(-1,1,10);
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
