function [vars, opts] = RPCA_R8_par1()
% Syntax:   [vars, opts] = RPCA_R8_par1();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.rawpath  = 'RPCA_R8_par1/data.mat';
vars.outpath  = 'RPCA_R8_par1.mat';
vars.r        = nan;
vars.lambdaL  = 0.525 * logspace(-1,1,11);
vars.lambdaS  = 0.010 * logspace(-1,1,11);
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
