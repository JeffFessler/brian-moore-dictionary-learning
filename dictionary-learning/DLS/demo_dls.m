%%
% Dictionary Least-Squares for photometric stereo
%
% Brian Moore
% brimoor@umich.edu
%

rng(1);

% Knobs
inpath = 'cat.mat';                         % Input images
SNR    = 30;                                % Image SNR, in DB
p      = 0.3;                               % Missing data probability
type   = '2c';                              % G2S model

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
Y        = Y .* Md;                         % Keep original mask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [I2N] IMAGES TO NORMALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NRMSE function
NRMSEI2N = @(X,Xtrue) norm(X(M3) - Xtrue(M3)) / norm(Xtrue(M3));

% True normals
Ntrue    = pinv(L) * I;                     % (3 x mn)

% LS
Nls      = pinv(L) * Y;                     % (3 x mn)
Nls_err  = NRMSEI2N(Nls,Ntrue);

% DLS
lambdaI2N        = 2.5;
muI2N            = 0.063;
optsI2N.A        = L;
optsI2N.xdim     = [3, m, n];
optsI2N.pdim     = [3, 8, 8];
optsI2N.pgap     = [1, 4, 4];
optsI2N.type     = 'soft';
optsI2N.nIters   = 10;
optsI2N.nItersDB = 1;
optsI2N.nItersX  = 25;
optsI2N.X0       = Nls;
optsI2N.Xtrue    = Ntrue;
optsI2N.NRMSEfcn = NRMSEI2N;
optsI2N.accel    = true;
optsI2N.flag     = 1;
[Ndls, ~, ~, statsI2N] = dls(Y,lambdaI2N,muI2N,optsI2N);
Ndls_err = statsI2N.nrmse(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [G2S] GRADIENTS TO SURFACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NRMSE function
NRMSEG2S = @(f,ftrue) norm(f(M) - ftrue(M)) / norm(ftrue(M));

% Construct G2S models
[~, ~    , ~    , ftrue] = G2Smodel(Ntrue,M,type);
[~, ~    , ~    , fls  ] = G2Smodel(Nls  ,M,type);
[A, bdls , normA, f0   ] = G2Smodel(Ndls ,M,type);
fls_err = NRMSEG2S(fls,ftrue);

% DLS-DLS
lambdaG2S        = 6.3;
muG2S            = 0.063;
optsG2S.A        = A;
optsG2S.xdim     = [m, n];
optsG2S.pdim     = [8, 8];
optsG2S.pgap     = [2, 2];
optsG2S.type     = 'soft';
optsG2S.nIters   = 10;
optsG2S.nItersDB = 1;
optsG2S.nItersX  = 25;
optsG2S.X0       = f0;
optsG2S.Xtrue    = ftrue;
optsG2S.NRMSEfcn = NRMSEG2S;
optsG2S.accel    = true;
optsG2S.tau      = (0.99 + ~optsG2S.accel) / normA^2;
optsG2S.flag     = 1;
[fdls, ~, ~, statsG2S] = dls(bdls,lambdaG2S,muG2S,optsG2S);
fdls_err = NRMSEG2S(fdls,ftrue);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape into matrices
Ftrue = reshape(ftrue,[m, n]) .* M;         % (m x n)
Fls   = reshape(fls  ,[m, n]) .* M;         % (m x n)
Fdls  = reshape(fdls ,[m, n]) .* M;         % (m x n)

%--------------------------------------------------------------------------
% Plot normals
%--------------------------------------------------------------------------
figure;

% True normals
subplot(1,3,1);
imshow(tilePatches(Ntrue',[m n],[3 1]),[]);
title('Ntrue');

% LS normals
subplot(1,3,2);
imshow(tilePatches(Nls',[m n],[3 1]),[]);
title(sprintf('Nls [NRMSE = %.2f]',Nls_err));

% DLS normals
subplot(1,3,3);
imshow(tilePatches(Ndls',[m n],[3 1]),[]);
title(sprintf('Ndls [NRMSE = %.2f]',Ndls_err));
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot surfaces
%--------------------------------------------------------------------------
figure;

% True surface
subplot(1,3,1);
surfplot(Ftrue);
title('Ftrue');

% LS surface
subplot(1,3,2);
surfplot(Fls);
title(sprintf('Fls [NRMSE = %.2f]',fls_err));

% Dictionary-LS surface
subplot(1,3,3);
surfplot(Fdls);
title(sprintf('Fdls [NRMSE = %.2f]',fdls_err));
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot optimization outputs
%--------------------------------------------------------------------------
figure;

% Cost function values
subplot(2,2,1);
plot(1:statsI2N.nIters,statsI2N.cost,'b-o');
xlabel('Iteration (k)');
ylabel('\Psi(N_k)');
title('Cost [I2N]');
axis tight; padAxis();

% NRMSE
subplot(2,2,2);
phndl    = zeros(1,2);
phndl(1) = plot(1:statsI2N.nIters,statsI2N.nrmse,'b-o'); hold on;
phndl(2) = plot([1, statsI2N.nIters],Nls_err  * [1, 1],'r-');
legend(phndl,'dLS','LS');
xlabel('Iteration (k)');
ylabel('NRMSE(N_k)');
title('NRMSE [I2N]');
axis tight; padAxis();

% Cost function values
subplot(2,2,3);
plot(1:statsG2S.nIters,statsG2S.cost,'b-o');
xlabel('Iteration (k)');
ylabel('\Psi(F_k)');
title('Cost [G2S]');
axis tight; padAxis();

% NRMSE
subplot(2,2,4);
phndl    = zeros(1,2);
phndl(1) = plot(1:statsG2S.nIters,statsG2S.nrmse,'b-o'); hold on;
phndl(2) = plot([1, statsG2S.nIters],fls_err  * [1 1],'r-');
legend(phndl,'dLS','LS');
xlabel('Iteration (k)');
ylabel('NRMSE(F_k)');
title('NRMSE [G2S]');
axis tight; padAxis();
%--------------------------------------------------------------------------
