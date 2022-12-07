function [vars, opts] = pincat_RPCA_R8_par1()
% Syntax:   [vars, opts] = pincat_RPCA_R8_par1();

% Knobs
vars.inpath   = 'pincat_R8.mat';
vars.rawpath  = 'pincat_RPCA_R8_par1/data.mat';
vars.outpath  = 'pincat_RPCA_R8_par1.mat';
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
