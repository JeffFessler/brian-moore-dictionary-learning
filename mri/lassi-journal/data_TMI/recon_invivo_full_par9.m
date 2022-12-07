function [vars, opts, vars0, opts0] = recon_invivo_full_par9()
% Syntax:   [vars, opts, vars0, opts0] = recon_invivo_full_par9();

% Knobs
vars.inpath   = 'invivo_full.mat';
vars.outpath  = 'recon_invivo_full_par9.mat';
vars.nLines   = 25; % [5, 10, 15, 20, 25, 30];
vars.seed     = 5;

% NO LASSI
opts = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ktSLR initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.p        = 1;
vars0.mu1      = 0;
vars0.mu2      = 0.001;
vars0.nItersO  = 5;
vars0.nItersI  = 10;
vars0.nItersCG = 5;

% ktSLR opts
%opts0.nItersO  = vars0.nItersO;
%opts0.nItersI  = vars0.nItersI;
%opts0.nItersCG = vars0.nItersCG;
opts0.beta1     = 1;
opts0.beta2     = 1;
opts0.beta1rate = 50;
opts0.beta2rate = 25;
opts0.stepSize  = [1, 1, 0.303];
