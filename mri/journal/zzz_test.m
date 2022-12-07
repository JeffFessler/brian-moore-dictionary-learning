%% L + S tests

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
lambdaL = 0.01;
lambdaS = 0.01;

% Model
A = 1;
normA = norm(A);

% Low-rank matrix
r = length(theta);
[U,~] = qr(randn(m,r),0);
[V,~] = qr(randn(n,r),0);
L = U * diag(theta) * V';

% Sparse matrix
S = sigmaS * sign(randn(m,n)) .* (rand(m,n) < p);

% Generate measurements
Y = A * (L + S) + (sigmaX / sqrt(max(m,n))) * randn(m,n);
%C = max(svd(A' * Y));
C = sqrt(max(m,n));

%--------------------------------------------------------------------------
% Proximal gradient method
%--------------------------------------------------------------------------
% Model
model = struct();
model.A = @(X) A * X;
model.At = @(X) A' * X;
model.lambdaL = C * lambdaL;
model.lambdaS = lambdaS;
model.r = r;

% Optimization parameters
opt = struct();
opt.tau = 0.95 / abs(normA)^2; % in [0 1] / norm(A)^2
opt.tol = 1e-6;
opt.Nmin = 5;
opt.Nmax = 100;
opt.method = method;

% Perform reconstruction
[Lhat_PGM Shat_PGM] = LpS_MRI(Y,model,opt);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ADMM method  
%--------------------------------------------------------------------------
% Model
model = struct();
model.A = @(X) A * X;
model.At = @(X) A' * X;
switch lower(method)
    case 'svt'
        params = {C * lambdaL};
    case 'optshrink'
        params = {r};
end
Lreg = struct('prox',method,'params',{params});
Sreg = struct('prox','soft','params',{{lambdaS}});
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
    minimize(0.5 * square_pos(norm(Y - A * (Lhat_CVX + Shat_CVX),'fro')) + (C * lambdaL) * norm_nuc(Lhat_CVX) + lambdaS * norm(vec(Shat_CVX),1));
cvx_end
%--------------------------------------------------------------------------

% Display results
NRMSE = @(X,Y) norm(X(:) - Y(:)) / norm(Y(:));
figure;
subplot(2,4,1); imagesc(L); colorbar; title('Ltrue');
subplot(2,4,2); imagesc(Lhat_CVX); colorbar; title(sprintf('Lhat (CVX) - [NRMSE = %.3f]',NRMSE(Lhat_CVX,L)));
subplot(2,4,3); imagesc(Lhat_PGM); colorbar; title(sprintf('Lhat (PGM) - [NRMSE = %.3f]',NRMSE(Lhat_PGM,L)));
subplot(2,4,4); imagesc(Lhat_ADMM); colorbar; title(sprintf('Lhat (ADMM) - [NRMSE = %.3f]',NRMSE(Lhat_ADMM,L)));
subplot(2,4,5); imagesc(S); colorbar; title('Strue');
subplot(2,4,6); imagesc(Shat_CVX); colorbar; title(sprintf('Shat (CVX) - [NRMSE = %.3f]',NRMSE(Shat_CVX,S)));
subplot(2,4,7); imagesc(Shat_PGM); colorbar; title(sprintf('Shat (PGM) - [NRMSE = %.3f]',NRMSE(Shat_PGM,S)));
subplot(2,4,8); imagesc(Shat_ADMM); colorbar; title(sprintf('Shat (ADMM) - [NRMSE = %.3f]',NRMSE(Shat_ADMM,S)));

%% L & S tests

rng(1);

% Generative parameters
m = 40;
n = 30;
theta = [5 2];
p = 0.1;
sigma = 0.1;

% Optimization parameters
method = 'SVT'; % {'SVT','OptShrink'}
lambdaL = 0.01;
lambdaS = 0.01;

% Model
A = 1;
normA = norm(A);

% Low-rank and sparse matrix
r = length(theta);
U = randn(m,r) .* (rand(m,r) < sqrt(p));
U = bsxfun(@times,1 ./ sqrt(sum(abs(U).^2,1)),U);
V = randn(n,r) .* (rand(n,r) < sqrt(p));
V = bsxfun(@times,1 ./ sqrt(sum(abs(V).^2,1)),V);
Xtrue = U * diag(theta) * V';

% Generate measurements
Y = A * Xtrue + (sigma / sqrt(max(m,n))) * randn(m,n);
%C = max(svd(A' * Y));
C = sqrt(max(m,n));

%{
%--------------------------------------------------------------------------
% ALM method
%--------------------------------------------------------------------------
% Construct algo model
model = struct();
model.A = @(X) A * X;
model.At = @(X) A' * X;
model.lambdaS = lambdaS;
model.X0 = A' * Y;

% Construct optimization parameters
opt = struct();
opt.method = method;
opt.beta0 = 1e6;
opt.betaInc = 25;
opt.Ninner = 10;
opt.Nouter = 10;

% Construct GD parameters
GDopt = struct();
GDopt.k = 1;
GDopt.normA = normA;
GDopt.tau = 1;

% Process based on L-method
switch lower(method)
    case 'svt'
        model.lambdaL = C * lambdaL;
    case 'optshrink'
        model.r = r;
end

% Run algorithm
Xhat_ALM = LaS_MRI(Y,model,opt,GDopt);
%--------------------------------------------------------------------------
%}

%--------------------------------------------------------------------------
% ADMM method  
%--------------------------------------------------------------------------
% Model
model = struct();
model.A = @(X) A * X;
model.At = @(X) A' * X;
switch lower(method)
    case 'svt'
        params = {C * lambdaL};
    case 'optshrink'
        params = {r};
end
Lreg = struct('prox',method,'params',{params});
Sreg = struct('prox','soft','params',{{lambdaS}});
model.X = {{Lreg Sreg}};
model.X0 = A' * Y;

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
Xhat_ADMM = X_ADMM{1};
%--------------------------------------------------------------------------

%{
%--------------------------------------------------------------------------
% CVX
%--------------------------------------------------------------------------
cvx_begin
    variables Xhat_CVX(m,n);
    minimize(0.5 * square_pos(norm(Y - A * (Xhat_CVX),'fro')) + (C * lambdaL) * norm_nuc(Xhat_CVX) + lambdaS * norm(vec(Xhat_CVX),1));
cvx_end
%--------------------------------------------------------------------------
%}

% Display results
NRMSE = @(X,Y) norm(X(:) - Y(:)) / norm(Y(:));
figure;
subplot(2,2,1); imagesc(Xtrue); colorbar; title('Xtrue');
%subplot(2,2,2); imagesc(Xhat_CVX); colorbar; title(sprintf('Xhat (CVX) - [NRMSE = %.3f]',NRMSE(Xhat_CVX,Xtrue)));
%subplot(2,2,3); imagesc(Xhat_ALM); colorbar; title(sprintf('Xhat (ALM) - [NRMSE = %.3f]',NRMSE(Xhat_ALM,Xtrue)));
subplot(2,2,4); imagesc(Xhat_ADMM); colorbar; title(sprintf('Xhat (ADMM) - [NRMSE = %.3f]',NRMSE(Xhat_ADMM,Xtrue)));

%% GenerateMovie tests

% Knobs
watermark = {'© Brian Moore - 2014';
             'University of Michigan'};
fps = 30;
res = 512;

% Generate test data
n = 128;
nt = 256;
x = (1:n)';
X = zeros(n,n,nt);
for i = 1:nt
    xi = circshift(x,i - 1);
    X(:,:,i) = xi * xi';
end

% Generate movie
%GenerateMovie(X,'test',fps,res);
GenerateMovie1(X,'test1',fps,watermark,res);
%GenerateMovie2(X,'test2',fps,watermark,res);

%% PlayMovie tests

% Knobs
path = 'xylophone.mpg'; % Path to movie file to load

% Read video frames
t = tic;
vidObj = VideoReader(path);
frames = double(read(vidObj));
fps = vidObj.FrameRate;
fprintf('Video loaded - (Time = %.2fs)\n',toc(t));

% Convert to grayscale frame "cube"
t = tic;
rgb_conv = permute([0.2989 0.5870 0.1140],[1 3 2 4]); % Conversion weights
X = squeeze(sum(bsxfun(@times,frames,rgb_conv),3));
[ny nx nt] = size(X);
fprintf('Grayscale conversion complete - (Time = %.2fs)\n',toc(t));

% Play movie
PlayMovie(X,fps);

%% Noise SNR test

% Knobs
SNR = 40;

Xtrue = truth.Xtrue(:,:,1);

SNR2sigma = @(SNR,X) exp(-SNR / 20) * norm(X(:)) / sqrt(numel(X)) / sqrt(2);
sigma = SNR2sigma(SNR,Xtrue) / sqrt(2);

[m n] = size(Xtrue);
Y = abs(Xtrue + sigma * (randn(m,n) + 1i * randn(m,n)))';

figure;
subplot(1,2,1);
imshow(abs(Xtrue)',[]);
subplot(1,2,2);
imshow(Y,[]);

%% svds w/ function handle test

m = 1000;
n = 500;

A = randn(m,n);

% Method #1: Direct
s1 = svd(A); s1 = s1(1) %#ok

% Method #2: svds() with function handle
AAtfcn = @(x) A' * (A * x);
opts = struct('issym',true,'isreal',false,'maxit',10,'disp',1);
s1c = sqrt(eigs(AAtfcn,n,1,'LM',opts)) %#ok

%% Generate Fessler data

% Knobs
SNR = 50;
Ncoil = 8;
Nframe = 20;
n_tr_merge = 100;
n_tr_round = 24;

% Generate data
ir_mri_dce_sim1;
Y = yi(:);
nd = nnz(mask(:));
nt = Nframe;
Xtrue = masker(xtrue,mask);

% Load truth (Xtrue480)
data = load('dynobj.mat');
Xtrue480 = data.dyn_obj;

% Non-iterative reconstruction
kdown = reshape(permute(yi,[1 3 2]),[],Ncoil);
kfull = data_share_fill(kdown,samp1);
xc = permute(ir_fftshift2(ifft2(ir_fftshift2(kfull))),[1 2 4 3]);
Xfft = squeeze(sum(xc .* conj(repmat(smap,[1 1 1 Nframe])),3)) ./ ...
       repmat(sum(abs(smap).^2,3),[1 1 Nframe]);
X0 = masker(Xfft,mask);

% NRMSE function
NRMSEfcn = @(Xhat) norm(Xhat(:) - Xtrue480(:)) / norm(Xtrue480(:));

% Spline interpolator
tt = ((1:Nframe) - 0.5) / Nframe * dce.duration_s / 60; % [min]
splineInterp = @(X) spline(tt,X,dce.ti);

%{
% Plot contrast curves
figure;
p = plot(dce.ti,dce.sig(labs,:));
xlabel('t');
ylabel('Brightness');
title('Contrast Curves');
legend(p,'ROI #1','ROI #2','ROI #3');
%}

%{
% Plot ground truth slices
figure;
subplot(1,3,1); imshow(xtrue(:,:,1)',[min(xtrue(:)) max(xtrue(:))]); title('Beginning');
subplot(1,3,2); imshow(xtrue(:,:,4)',[min(xtrue(:)) max(xtrue(:))]); title('Peak Intensity');
subplot(1,3,3); imshow(xtrue(:,:,20)',[min(xtrue(:)) max(xtrue(:))]); title('End');
%}

%% Perform L + S

% Knobs
algoStr = 'L+S PGM'; % {'L+S PGM','L+S ADMM'}
methodStr = 'OptShrink'; % {'SVT','OptShrink'}
lambdaL = 0.01; % pct of /sigma_{r+1}        [0.001 1]
lambdaS = 0.01; % pct of \|x_{true}\|_{\inf} [0.001 0.05]
r = 1;
%nb = [24 27];

% System matrix norm
%AAtfcn = @(x) A' * (A * x);
%opts = struct('issym',true,'isreal',false,'maxit',5,'disp',1);
%normA = sqrt(eigs(AAtfcn,size(A,2),1,'LM',opts)) %#ok
normA = 20.526;

% Data
T = 1; % TempFFT(2);
C = svd(Xtrue); C = C(r + 1);
Xtrue_inf = max(abs(Xtrue(:)));

% Construct algo arguments
method.normA = normA;
method.str = methodStr;
switch lower(methodStr)
    case 'svt'
        method.lambdaL = C * lambdaL;
        if exist('nb','var')
            Nb = floor(nx / nb(1)) * floor(ny / nb(2));
            method.lambdaL = model.lambdaL / Nb;
        end
    case 'optshrink'
        method.r = r;
end
method.lambdaS = Xtrue_inf * lambdaS;
if exist('nb','var')
    method.nb = nb;
end
method.mask = mask;
method.X0 = X0;

% Generate true data
truth.Xtrue = Xtrue;

% Run algo
[NRMSE recon] = RunLSAlgo(algoStr,Y,A,T,nd,1,nt,method,truth);
Lhat = embed(recon.L,mask);
Shat = embed(recon.S,mask);
Xhat = embed(recon.X,mask);

% Interpolate results
Xfft480 = splineInterp(abs(Xfft));
Xhat480 = splineInterp(abs(Xhat));
Lhat480 = splineInterp(abs(Lhat));
Shat480 = splineInterp(abs(Shat));

% NRMSEs
NRMSE_fft_LpS = [NRMSEfcn(Xfft480) NRMSEfcn(Xhat480)] %#ok

%% Display L + S results

% Concatenate slabs
M = {Xtrue480 Xfft480 Xhat480 Lhat480 Shat480};  
for i = 1:numel(M)
    % Permute rows/columns
    M{i} = permute(M{i},[2 1 3]);
    
    % Display each on full-scale
    %M{i} = M{i} - min(M{i}(:));
    %M{i} = M{i} / max(M{i}(:));
    
    % Display on Xtrue scale
    M{i} = M{i} - min(Xtrue480(:));
    M{i} = M{i} / max(Xtrue480(:));
    M{i} = min(max(M{i},0),1);
end
M = cat(2,M{:});

% Play movie
c = 1.5;
fps = 30;
PlayMovie(M,fps,[nan round(c * size(M,1))]);

%% Perform L & S

% Knobs
algoStr = 'L&S ADMM';
methodStr = 'OptShrink'; % {'SVT','OptShrink'}
lambdaL = 0.5; % pct of /sigma_{r+1}        [0.001 1]
lambdaS = 0.005; % pct of \|x_{true}\|_{\inf} [0.001 0.05]
r = 4;
%nb = [24 27]; % {8,24}

% Dimensions
nd = nnz(mask(:));
nt = Nframe;

% System matrix norm
%AAtfcn = @(x) A' * (A * x);
%opts = struct('issym',true,'isreal',false,'maxit',5,'disp',1);
%normA = sqrt(eigs(AAtfcn,nd * nt,1,'LM',opts)) %#ok
normA = 20.526;

% Data
T = TempFFT(2);
C = svd(Xtrue); C = C(r + 1);
Xtrue_inf = max(abs(Xtrue(:)));

% Construct algo arguments
method.normA = normA;
method.str = methodStr;
switch lower(methodStr)
    case 'svt'
        method.lambdaL = C * lambdaL;
        if exist('nb','var')
            Nb = floor(nx / nb(1)) * floor(ny / nb(2));
            method.lambdaL = model.lambdaL / Nb;
        end
    case 'optshrink'
        method.r = r;
end
method.lambdaS = Xtrue_inf * lambdaS;
if exist('nb','var')
    method.nb = nb;
end
method.X0 = X0;

% Generate true data
truth.Xtrue = Xtrue;

% Run algo
[NRMSE recon] = RunLSAlgo(algoStr,Y,A,T,nd,1,nt,method,truth);
Xhat = embed(recon.X,mask);

% Interpolate results
Xfft480 = splineInterp(abs(Xfft));
Xhat480 = splineInterp(abs(Xhat));

% NRMSEs
NRMSE_fft_LaS = [NRMSEfcn(Xfft480) NRMSEfcn(Xhat480)] %#ok

%% Display L & S results

% Concatenate slabs
M = {Xtrue480 Xfft480 Xhat480};
for i = 1:numel(M)
    % Permute rows/columns
    M{i} = permute(M{i},[2 1 3]);
    
    % Display each on full-scale
    %M{i} = M{i} - min(M{i}(:));
    %M{i} = M{i} / max(M{i}(:));
    
    % Display on Xtrue scale
    M{i} = M{i} - min(Xtrue480(:));
    M{i} = M{i} / max(Xtrue480(:));
    M{i} = min(max(M{i},0),1);
end
M = cat(2,M{:});

% Play movie
c = 1.5;
fps = 30;
PlayMovie(M,fps,[NaN round(c * size(M,1))]);

%% Raj Meeting 2/20/15
%
% Check to see if lesions are identifiable in FFT reconstruction
% 

% First run "Generate Fessler data" cell

Xfft480 = splineInterp(abs(Xfft));
Mm = reshape(Xfft480,[],480);

[Um Sm Vm] = svd(Mm,'econ');
M0 = Sm(1) * Um(:,1) * Vm(:,1)';

XX = Mm - M0;

figure;

subplot(1,2,1);
v1 = sum(Mm - M0,2);
semilogy(v1);

subplot(1,2,2);
v2 = sum((Mm - M0).^2,2);
semilogy(v2);


[~,idx] = sort(v2,'descend');

[i j] = ind2sub([192 216],idx(1:3));

figure;
plot(j,i,'bo');

u1 = squeeze(Xfft480(i(1),j(1),:));
u2 = squeeze(Xfft480(i(2),j(2),:));
u3 = squeeze(Xfft480(i(3),j(3),:));

figure;
subplot(1,3,1);
plot(u1,'r-o');
subplot(1,3,2);
plot(u2,'r-o');
subplot(1,3,3);
plot(u3,'r-o');
