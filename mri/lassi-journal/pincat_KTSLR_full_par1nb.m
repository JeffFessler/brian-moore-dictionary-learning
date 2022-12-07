function [vars, opts] = pincat_KTSLR_full_par1nb()
% Syntax:   [vars, opts] = pincat_KTSLR_full_par1nb();

% Knobs
vars.inpath   = 'pincat_full.mat';
vars.rawpath  = 'pincat_KTSLR_full_par1nb/data.mat';
vars.outpath  = 'pincat_KTSLR_full_par1nb.mat';
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = 1:5;
vars.SNR      = 45; % Add noise
vars.p        = 1;
vars.mu1      = 2;
vars.mu2      = 0.01;
vars.nItersO  = 5;
vars.nItersI  = 10;
vars.nItersCG = 5;

% ktSLR opts
%opts.nItersO  = vars.nItersO;
%opts.nItersI  = vars.nItersI;
%opts.nItersCG = vars.nItersCG;
opts.beta1     = 1;
opts.beta2     = 1;
opts.beta1rate = 50;
opts.beta2rate = 25;
opts.stepSize  = [1, 1, 0.303];
%opts.Xtrue    = data.Xtrue;
