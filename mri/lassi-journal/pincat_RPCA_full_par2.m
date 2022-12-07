function [vars, opts] = pincat_RPCA_full_par2()
% Syntax:   [vars, opts] = pincat_RPCA_full_par2();

% Knobs
vars.inpath   = 'pincat_full.mat';
vars.rawpath  = 'pincat_RPCA_full_par2/data.mat';
vars.outpath  = 'pincat_RPCA_full_par2.mat';
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = 1:5;
vars.r        = 9;
vars.lambdaL  = nan;
vars.lambdaS  = 0.005;
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
