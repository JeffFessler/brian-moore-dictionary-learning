function [vars, opts, vars0, opts0] = recon_otazo_full_par12()
% Syntax:   [vars, opts, vars0, opts0] = recon_otazo_full_par12();

% Knobs
vars.inpath    = 'otazo_full.mat';
vars.outpath   = 'recon_otazo_full_par12.mat';
vars.ps        = 1 ./ 20; % [4, 8, 12, 16, 20, 24]
vars.seed      = 3;

% NO LASSI
opts = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ktSLR initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vars0.p        = 1;
vars0.mu1      = 0.2;
vars0.mu2      = 0.003;
vars0.nItersO  = 5;
vars0.nItersI  = 10;
vars0.nItersCG = 5;

% ktSLR opts
%opts0.nItersO  = vars.nItersO;
%opts0.nItersI  = vars.nItersI;
%opts0.nItersCG = vars.nItersCG;
opts0.beta1     = 1;
opts0.beta2     = 1;
opts0.beta1rate = 50;
opts0.beta2rate = 25;
opts0.stepSize  = [1, 1, 0.303];
%opts0.Xtrue    = data.Xtrue;
