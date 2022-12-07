function [vars, opts] = pincat_KTSLR_R8_par1n()
% Syntax:   [vars, opts] = pincat_KTSLR_R8_par1n();

% Knobs
vars.inpath   = 'pincat_R8n.mat';
vars.rawpath  = 'pincat_KTSLR_R8_par1n/data.mat';
vars.outpath  = 'pincat_KTSLR_R8_par1n.mat';
vars.p        = [0.1, 0.5, 1];
vars.mu1      = 0.0200 * logspace(-2,2,9);
vars.mu2      = 0.0001 * logspace(-2,2,9);
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
