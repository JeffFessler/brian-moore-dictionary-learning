function [vars, opts] = pincat_RPCA_R8_par2n()
% Syntax:   [vars, opts] = pincat_RPCA_R8_par2n();

% Knobs
vars.inpath   = 'pincat_R8n.mat';
vars.rawpath  = 'pincat_RPCA_R8_par2n/data.mat';
vars.outpath  = 'pincat_RPCA_R8_par2n.mat';
vars.r        = 1:20;
vars.lambdaL  = nan;
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
