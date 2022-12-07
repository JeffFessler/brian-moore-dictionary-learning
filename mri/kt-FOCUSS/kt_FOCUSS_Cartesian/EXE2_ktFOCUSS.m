%% k-t FOCUSS

% Knobs
inpath = 'Y.mat';
factor = 0.5;
lambda = 0;
Mouter = 2;
Minner = 40;

% Load data
load(inpath);
Y    = DownSino;
mask = double(mask);
[nx, ny, nt] = size(Y);

% Discard high frequency components
Low_Y = Y;
Low_Y((0.5 * num_low_phase + 1):(nx - 0.5 * num_low_phase),:,:) = 0;

% Add dependencies
addpath('bin');
addpath('data');

% System matrix
A  = @(x,mask) fft(x,[],1) .* mask;
AT = @(x,mask) ifft(x .* mask,[],1);

% k-t FOCUSS
Xhat = KTFOCUSS(A,AT,Y,Low_Y,mask,factor,lambda,Minner,Mouter);

% Display results
PlayMovie(Xhat);
