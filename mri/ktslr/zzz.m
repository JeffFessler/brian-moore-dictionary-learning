%% Check data generation

% Knobs
nLines  = 20;
SNR     = inf;
seed    = 1;
%inpath = 'invivo_full.mat';
inpath  = 'pincat_full.mat';

% Generate invivo data
[Y, ~, ~, Xtrue, X0] = generateInvivoData(nLines,SNR,seed,inpath);
mask = (abs(Y) ~= 0);

% Visualize data
PlayMovie(cat(2,mask,Xtrue,X0));

%% Generate ~R8 data

% Knobs
%{
nLines  = 15; % invivo
SNR     = inf;
seed    = 1;
inpath  = 'invivo_full.mat';
outpath = 'invivo_R8.mat';
%}
%{
nLines   = 18;
SNR      = inf;
seed     = 1;
inpath   = 'pincat_full.mat';
outpath  = 'pincat_R8.mat';
%}
nLines   = 18;
SNR      = 45;
seed     = 1;
inpath   = 'pincat_full.mat';
outpath  = 'pincat_R8n.mat';

% Generate invivo data
[Y, ~, ~, Xtrue, Xfft] = generateInvivoData(nLines,SNR,seed,inpath);
mask = (abs(Y) ~= 0);

% View mask
ps = nnz(mask) / numel(mask); % Empirical sampling factor
opts.xlabels = sprintf('ps = %.3f',ps);
PlayMovie(fftshift(mask),opts);

% Save undersampled data
if exist('outpath','var') && ~isempty(outpath)
    save(outpath,'Xfft','Xtrue','Y','mask');
    fprintf('File "%s" written\n',outpath);
end

%% Validate Afft

% Knobs
[ny, nx, nt] = deal(128,128,50);
nLines = 15;

% Sampling mask
mask = strucrand(ny,nx,nt,nLines);
mask = logical(fftshift(fftshift(mask,1),2));

% Afft
siz = [ny, nx, nt];
A0  = Afft(mask,siz);

% ktslr
S   = find(mask);
A1  = @(z)  A_fhp3D(z,S);
A1t = @(z) reshape(At_fhp3D(z,S,ny,nx,nt),[],nt);

% Test forward operation
X = randn(siz);
A0X = A0 * X;
A1X = A1(X);
sizeA0X = size(A0X) %#ok
sizeA1X = size(A1X) %#ok
errA = norm(A0X(:) - A1X(:)) %#ok

% Test backward operation
A0tA0X = A0' * A0X;
A1tA1X = A1t(A1X);
sizeA0tA0X = size(A0tA0X) %#ok
sizeA1tA1X = size(A1tA1X) %#ok
errAt = norm(A0tA0X(:) - A1tA1X(:)) %#ok
