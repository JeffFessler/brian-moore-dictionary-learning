%%
% Demo various G2S formulations
%
% Brian Moore
% brimoor@umich.edu
%

rng(1);

% Data knobs
inpath = 'cat.mat';     % Input images

% Fast solver knobs
typef  = '2a';          % Fast solver type

% Iterative solver knobs
typei  = '1n';          % Iterative solver type
nIters = 1000;          % # iterations
accel  = true;          % Accelerated LS?

% Load data
load(inpath);

% Normals
[m, n, d] = size(I);
N = reshape(reshape(I,m * n,d) * pinv(L),[m, n, 3]);

% Fast solver
[~, ~, ~, ffast] = G2Smodel(N,M,typef);
Ffast            = reshape(ffast,m,n);

% Iterative solver
[A, b, normA] = G2Smodel(N,M,typei);
opts.nIters   = nIters;
opts.accel    = accel;
opts.tau      = (0.99 + ~accel) / normA^2;
opts.flag     = 0;
fils          = ils(A,b,opts);

% Optimize offset
c    = fminsearch(@(c) norm(ffast - (fils - c)),0) %#ok
fils = fils - c;
Fils = reshape(fils,[m, n]);

% Display error
err = norm(ffast - fils) %#ok

% Display results
figure;
subplot(1,2,1); surfplot(Ffast); title('Ffast');
subplot(1,2,2); surfplot(Fils);  title('Fils');
