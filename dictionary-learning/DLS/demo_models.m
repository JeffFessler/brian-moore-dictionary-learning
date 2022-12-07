%%
% Compare gradients-to-surface (G2S) models
%
% Brian Moore
% brimoor@umich.edu
%

rng(1);

% Knobs
inpath = 'cat.mat';                         % Input images
SNR    = 50;                                % Image SNR, in DB
p      = 0.05;                              % Missing data probability
zlim   = [-150, 150];                       % Clipped z-limits (plotting)

% Load data
load(inpath);

% Reshape data
[m, n, d] = size(I);
Md = repmat(M(:),1,d)';                     % (d x mn)
M3 = repmat(M(:),1,3)';                     % (3 x mn)
I  = reshape(I,m * n,d)';                   % (d x mn)
L  = L';                                    % (d x 3)

% Corrupt images
SIGMA    = exp(-SNR / 20.0) * norm(I(:)) / sqrt(nnz(Md));
ZEROS    = (rand(size(I)) < p);
Y        = I + SIGMA * randn(size(I));      % Add noise
Y(ZEROS) = 0;                               % Add missing data
Y        = Y .* Md;                         % Apply original mask

% Compute normals from corrupted images via least-squares
Nls      = pinv(L) * Y;                     % (3 x mn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform G2S
% All of these minimize \|b - Af\|_2^2 for different (A,b)
% See G2Smodel for more detail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order model with circulant boundary conditions
% Minimize via gradient descent
[A, b, normA] = G2Smodel(Nls,M,'1c');
opts.nIters = 1000;
opts.accel  = true;                         % Nesterov acceleration
opts.tau    = (0.99 + ~opts.accel) / normA^2;
opts.flag   = 0;
f1c = ils(A,b,opts);

% Second order "integrability" model
% Minimized in closed-form via 2D DST
[~, ~, ~, f2a] = G2Smodel(Nls,M,'2a');

% Second order "integrability" model
% Minimized in closed-form via 2D DCT
[~, ~, ~, f2c] = G2Smodel(Nls,M,'2c');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape into surfaces
pretty = @(f) max(zlim(1),min(zlim(2),reshape(f, [m, n])));
F1c = pretty(f1c);                          % (m x n)
F2a = pretty(f2a);                          % (m x n)
F2c = pretty(f2c);                          % (m x n)

% Plot surfaces
figure;
subplot(1,3,1); surfplot(F1c); title('F1c');
subplot(1,3,2); surfplot(F2a); title('F2a');
subplot(1,3,3); surfplot(F2c); title('F2c');
