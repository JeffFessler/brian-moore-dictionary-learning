%%
% Least-squares photometric stereo on noisy data
%
% Brian Moore
% brimoor@umich.edu
%

rng(42);

% Knobs
inpath = 'cat.mat';                         % Input images
SNR    = 1;                                 % SNR, in dB
type   = '2a';                              % Reconstruction method
zlim   = [-150, 150];                       % Plotting limits

% Load data
load(inpath);

% Reshape data
[m, n, d] = size(I);
Md = repmat(M(:),1,d)';                     % (d x mn)
M3 = repmat(M(:),1,3)';                     % (3 x mn)
I  = reshape(I,m * n,d)';                   % (d x mn)
L  = L';                                    % (d x 3)

% Add Poisson noise with desired SNR
alpha = (sum(I(:)) / sum(I(:).^2)) * 10^(0.1 * SNR);
Y     = poissrnd(alpha * I) / alpha;        % (d x mn)

%{
% Add Gaussian noise with desired SNR
sigma  = 10^(-SNR / 20) * (norm(I(Md)) / sqrt(nnz(Md)));
Y      = (I + sigma * randn(size(I))) .* Md;
%Y     = max(Y,0); % This helps, ALOT
%}

% Verify SNR
SNRact = 10 * log10(mean(I(:).^2) / mean((Y(:) - I(:)).^2)) %#ok

% Normals
Ntrue = pinv(L) * I;                        % (3 x mn)
Nhat  = pinv(L) * Y;                        % (3 x mn)

% Surfaces
[~, ~, ~, ftrue] = G2Smodel(Ntrue,M,type);  % (mn x 1)
[~, ~, ~, fhat]  = G2Smodel(Nhat,M,type);   % (mn x 1)
Ftrue = prettySurf(ftrue,M,zlim);           % (m x n)
Fhat  = prettySurf(fhat,M,zlim);            % (m x n)

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
idx = [1, 2];
figure('Color','w');

% Raw images
Id = zeros(m,0);
Yd = zeros(m,0);
for i = 1:numel(idx)
    Id = [Id, reshape(I(idx(i),:),m,n)]; %#ok
    Yd = [Yd, reshape(Y(idx(i),:),m,n)]; %#ok
end

% Plot results
subplot(2,2,1); imshow(Id,[0, 1]); title('I');
subplot(2,2,2); imshow(Yd,[0, 1]); title(sprintf('Y [%d dB]',SNR));
subplot(2,2,3); surfplot(Ftrue); box on; title('Ftrue');
subplot(2,2,4); surfplot(Fhat);  box on; title('Fhat');

% Save figure
%export_fig -png -m2 -nocrop pois
%export_fig -png -m2 -nocrop gaus_with_negs
%--------------------------------------------------------------------------
