function [vars, opts] = invivo_RPCA_full_par1()
% Syntax:   [vars, opts] = invivo_RPCA_full_par1();

% Knobs
vars.inpath   = 'invivo_full.mat';
vars.rawpath  = 'invivo_RPCA_full_par1/data.mat';
vars.outpath  = 'invivo_RPCA_full_par1.mat';
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = 1:5;
vars.r        = nan;
vars.lambdaL  = 0.1;
vars.lambdaS  = 0.001;
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
