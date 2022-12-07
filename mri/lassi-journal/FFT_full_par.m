function vars = FFT_full_par()
% Syntax:   vars = FFT_full_par();

% Knobs
vars.inpath    = 'otazo_full.mat';
vars.rawpath   = 'FFT_full_par/data.mat';
vars.outpath   = 'FFT_full_par.mat';
vars.ps        = 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed      = 1:5;
