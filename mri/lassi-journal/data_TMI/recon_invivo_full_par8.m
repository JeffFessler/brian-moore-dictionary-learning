function [vars, opts, vars0, opts0] = recon_invivo_full_par8()
% Syntax:   [vars, opts, vars0, opts0] = recon_invivo_full_par8();

% Knobs
vars.inpath   = 'invivo_full.mat';
vars.outpath  = 'recon_invivo_full_par8.mat';
vars.nLines   = 25; % [5, 10, 15, 20, 25, 30];
vars.seed     = 5;

% NO LASSI
opts = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPCA initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vars0.r        = nan;
vars0.lambdaL  = 0.158;
vars0.lambdaS  = 0.0025;
vars0.nIters   = 250;

% RPCA opts
%opts0.A      = A;
%opts0.T      = T;
%opts0.r      = r;
%opts0.nIters = vars.nIters;
%opts0.L0     = reshape(data.Xfft,[],nt);
%opts0.S0     = zeros(ny * nx,nt);
%opts0.Xtrue  = reshape(data.Xtrue,[],nt);
opts0.accel   = false;
opts0.tau     = 1;
opts0.flag    = 1;
