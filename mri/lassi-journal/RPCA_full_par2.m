function [vars, opts] = RPCA_full_par2()
% Syntax:   [vars, opts] = RPCA_full_par2();

% Knobs
vars.inpath   = 'otazo_full.mat';
vars.rawpath  = 'RPCA_full_par2/data.mat';
vars.outpath  = 'RPCA_full_par2.mat';
vars.ps       = 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed     = 1:5;
vars.r        = 1;
vars.lambdaL  = nan;
vars.lambdaS  = 0.01;
vars.nIters   = 250;

% RPCA opts
%opts.A       = A;
%opts.T       = T;
%opts.r       = r;
%opts.nIters  = vars.nIters;
%opts.L0      = reshape(data.Xfft,[],nt);
%opts.S0      = zeros(ny * nx,nt);
%opts.Xtrue   = reshape(data.Xtrue,[],nt);
opts.accel    = false;
opts.tau      = 1;
opts.flag     = 1;
