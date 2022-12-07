function [vars, opts] = pincat_RPCA_full_par2n()
% Syntax:   [vars, opts] = pincat_RPCA_full_par2n();

% Knobs
vars.inpath   = 'pincat_full.mat';
vars.rawpath  = 'pincat_RPCA_full_par2n/data.mat';
vars.outpath  = 'pincat_RPCA_full_par2n.mat';
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = 1:5;
vars.SNR      = 45; % Add noise
vars.r        = 1;
vars.lambdaL  = nan;
vars.lambdaS  = 0.0003;
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
