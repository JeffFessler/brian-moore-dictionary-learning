function [vars, opts] = RPCA_full_par1()
% Syntax:   [vars, opts] = RPCA_full_par1();

% Knobs
vars.inpath   = 'otazo_full.mat';
vars.rawpath  = 'RPCA_full_par1/data.mat';
vars.outpath  = 'RPCA_full_par1.mat';
vars.ps       = 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed     = 1:5;
vars.r        = nan;
vars.lambdaL  = 1.1955; % visual
%vars.lambdaL = 0.53;   % mmse
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
