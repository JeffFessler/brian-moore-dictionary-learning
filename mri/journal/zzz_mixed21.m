%% L + S with mixed 2-1 sparsity

rng(1);

% Knobs
rows = 1:200;
cols = 1:200;

% Load data
load('C:\Users\brimoor\Documents\Google Drive\MATLAB\Michigan Research\Summer 2013 - RMT Bulk Detection\OptShrink\data sets\yahoo movies\yahoo_movie_ratings_test.mat');
Y = M_test'; % n x m (w/ corrupted columns?)
clear M_test;

% Optimization parameters
method = 'OptShrink'; % {'SVT','OptShrink'}
lambdaL = 0.01;
lambdaS = 0.01;
r = 35;
normA = 1; % A == 1

lambdaL = max(svd(Y)) * lambdaL;
lambdaS = max(abs(Y(:))) * lambdaS;

%--------------------------------------------------------------------------
% Proximal gradient method
%--------------------------------------------------------------------------
% Model
model = struct();
model.A = @(X) X;
model.At = @(Y) Y;
model.lambdaL = lambdaL;
model.lambdaS = lambdaS;
model.r = r;

% Optimization parameters
opt = struct();
opt.tau = 0.95 / abs(normA)^2; % in [0 1] / norm(A)^2
opt.tol = 1e-6;
opt.Nmin = 5;
opt.Nmax = 20;
opt.method = method;

% Perform reconstruction
[Lhat_PGM Shat_PGM] = LpS21_MRI(Y,model,opt);
%[Lhat_PGM Shat_PGM] = LpS_MRI(Y,model,opt);
%--------------------------------------------------------------------------

% Display results
figure;
subplot(1,3,1); imagesc(Y(rows,cols)); colorbar; title('Y');
subplot(1,3,2); imagesc(Lhat_PGM(rows,cols)); colorbar; title('Lhat');
subplot(1,3,3); imagesc(Shat_PGM(rows,cols)); colorbar; title('Shat');

%% L + S with mixed 2-1 sparsity

rng(1);

% Generative parameters
m = 40;
n = 30;
theta = [5 2];
p = 0.1;
sigmaS = 1;
sigmaX = 0.1;

% Optimization parameters
method = 'SVT'; % {'SVT','OptShrink'}
lambdaL = 1;
lambdaS = 0.5;

% Model
A = 1;
normA = norm(A);

% Low-rank matrix
r = length(theta);
[U,~] = qr(randn(m,r),0);
[V,~] = qr(randn(n,r),0);
L = U * diag(theta) * V';

% Column-sparse matrix
S = sigmaS * (randn(m,n) .* repmat(rand(1,n) < p,[m 1]));

% Generate measurements
Y = A * (L + S) + (sigmaX / sqrt(max(m,n))) * randn(m,n);

%--------------------------------------------------------------------------
% Proximal gradient method
%--------------------------------------------------------------------------
% Model
model = struct();
model.A = @(X) A * X;
model.At = @(Y) A' * Y;
model.lambdaL = lambdaL;
model.lambdaS = lambdaS;
model.r = r;

% Optimization parameters
opt = struct();
opt.tau = 0.95 / abs(normA)^2; % in [0 1] / norm(A)^2
opt.tol = 1e-6;
opt.Nmin = 5;
opt.Nmax = 1000;
opt.method = method;

% Perform reconstruction
[Lhat_PGM Shat_PGM] = LpS21_MRI(Y,model,opt);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ADMM
%--------------------------------------------------------------------------
% Model
model = struct();
model.A = @(X) A * X;
model.At = @(Y) A' * Y;
switch lower(method)
    case 'svt'
        params = {lambdaL};
    case 'optshrink'
        params = {r};
end
Lreg = struct('prox',method,'params',{params});
Sreg = struct('prox','mixed21','params',{{lambdaS}});
model.X = {{Lreg} {Sreg}};

% Optimization parameters
opt = struct();
opt.rho = 1;
opt.epsilon = 1e-6;
opt.Nmin = 5;
opt.Nmax = 100;
opt.nu = 0.999;

% Construct GD parameters
GDopt = struct();
GDopt.k = 5;
GDopt.normA = normA;
GDopt.tau = 0.95; % in [0 2]

% Run algorithm
%X_ADMM = MSM_MRI(Y,model,opt,GDopt);
X_ADMM = MSMn_MRI(Y,model,opt,GDopt);
Lhat_ADMM = X_ADMM{1};
Shat_ADMM = X_ADMM{2};
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% CVX
%--------------------------------------------------------------------------
cvx_begin
    variables Lhat_CVX(m,n) Shat_CVX(m,n);
    minimize(0.5 * square_pos(norm(Y - A * (Lhat_CVX + Shat_CVX),'fro')) + lambdaL * norm_nuc(Lhat_CVX) + lambdaS * norm21(Shat_CVX));
cvx_end
%--------------------------------------------------------------------------

% Display results
NRMSE = @(X,Y) norm(X(:) - Y(:)) / norm(Y(:));
figure;
subplot(2,4,1); imagesc(L); colorbar; title('Ltrue');
subplot(2,4,2); imagesc(Lhat_PGM); colorbar; title(sprintf('Lhat (PGM) - [NRMSE = %.3f]',NRMSE(Lhat_PGM,L)));
subplot(2,4,3); imagesc(Lhat_ADMM); colorbar; title(sprintf('Lhat (ADMM) - [NRMSE = %.3f]',NRMSE(Lhat_ADMM,L)));
subplot(2,4,4); imagesc(Lhat_CVX); colorbar; title(sprintf('Lhat (CVX) - [NRMSE = %.3f]',NRMSE(Lhat_CVX,L)));
subplot(2,4,5); imagesc(S); colorbar; title('Strue');
subplot(2,4,6); imagesc(Shat_PGM); colorbar; title(sprintf('Shat (PGM) - [NRMSE = %.3f]',NRMSE(Shat_PGM,S)));
subplot(2,4,7); imagesc(Shat_ADMM); colorbar; title(sprintf('Shat (ADMM) - [NRMSE = %.3f]',NRMSE(Shat_ADMM,S)));
subplot(2,4,8); imagesc(Shat_CVX); colorbar; title(sprintf('Shat (CVX) - [NRMSE = %.3f]',NRMSE(Shat_CVX,S)));
