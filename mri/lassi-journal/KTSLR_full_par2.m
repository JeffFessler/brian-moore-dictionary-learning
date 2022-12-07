function [vars, opts] = KTSLR_full_par2()
% Syntax:   [vars, opts] = KTSLR_full_par2();

% Knobs
vars.inpath   = 'otazo_full.mat';
vars.rawpath  = 'KTSLR_full_par2/data.mat';
vars.outpath  = 'KTSLR_full_par2.mat';
vars.ps       = 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed     = 1:5;
vars.p        = 1;
vars.mu1      = 0.2;
vars.mu2      = 0.003;
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
