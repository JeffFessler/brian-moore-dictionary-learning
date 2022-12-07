function [vars, opts, vars0, opts0] = recon_pincat_full_par3()
% Syntax:   [vars, opts, vars0, opts0] = recon_pincat_full_par3();

% Knobs
vars.inpath   = 'pincat_full.mat';
vars.outpath  = 'recon_pincat_full_par3.mat';
vars.nLines   = 15; % [5, 10, 15, 20, 25, 30];
vars.seed     = 4;
vars.SNR      = 45; % Add noise

% No LASSI
opts = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPCA initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.r       = nan;
vars0.lambdaL = 3.793;
vars0.lambdaS = 0.15;
vars0.nIters  = 250;

% RPCA opts
%opts0.A      = A;
%opts0.T      = T;
%opts0.r      = vars0.r;
%opts0.nIters = vars0.nIters;
%opts0.L0     = reshape(Xfft,[],nt);
%opts0.S0     = zeros(ny * nx,nt);
opts0.accel   = false;
opts0.tau     = 1;
opts0.flag    = 1;
