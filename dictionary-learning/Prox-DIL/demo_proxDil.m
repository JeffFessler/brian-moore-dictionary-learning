%%
% Compare proxDil and soupDil
%
% Brian Moore
% brimoor@umich.edu
%

rng(42);

% Data knobs
d   = 64;
m   = 64;
n   = 4096;
p   = 0.10;
SNR = 40;

% Algorithm knobs
mu     = 0.3;           % Sparsity regularization
nIters = 100;           % # iterations
type   = 'hard';        % Type of thresholding
flag   = 0;             % Print iteration status?

% Generate data
Dtrue = randn(d,m);
Dtrue = bsxfun(@rdivide,Dtrue,sqrt(sum(abs(Dtrue).^2,1)));
Btrue = randn(m,n) .* (rand(m,n) < p);
Y = corrupt(Dtrue * Btrue,SNR);

% Initialization
D0 = dctmtx(max(d,m))'; D0 = D0(1:d,1:m);
B0 = zeros(m,n);

% SOUP-DIL
opts1 = struct();
opts1.D0     = D0;
opts1.B0     = B0;
opts1.nIters = nIters;
opts1.type   = type;
opts1.tau    = 0.0175; % Sparsity optimization
opts1.flag   = flag;
[D1, B1, stats1] = soupDil(Y,mu,opts1);

% Prox-DIL
opts2 = struct();
opts2.D0     = D0;
opts2.B0     = B0;
opts2.fixedD = false;
opts2.nIters = nIters;
opts2.type   = type;
opts2.accel  = true;
opts2.tau    = (0.99 + ~opts2.accel);
opts2.flag   = flag;
[D2, B2, stats2] = proxDil(Y,mu,opts2);

% Errors
NRMSEfcn = @(X,Xtrue) norm(X(:) - Xtrue(:)) / norm(Xtrue(:));
nrmseD1 = NRMSEfcn(D1,Dtrue) %#ok
nrmseB1 = NRMSEfcn(B1,Btrue) %#ok
nrmseD2 = NRMSEfcn(D2,Dtrue) %#ok
nrmseB2 = NRMSEfcn(B2,Btrue) %#ok

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure();
cm = linspecer(2);

% cost
subplot(2,2,1);
phndl = zeros(1,2);
phndl(1) = semilogy(1:nIters,stats1.cost,'-o','Color',cm(1,:)); hold on;
phndl(2) = semilogy(1:nIters,stats2.cost,'-o','Color',cm(2,:));
xlabel('iteration');
title('cost');
legend(phndl,'SOUP-DIL','Prox-DIL');
axis tight; padAxis();

% time
subplot(2,2,2);
plot(1:nIters,cumsum(stats1.time),'-o','Color',cm(1,:)); hold on;
plot(1:nIters,cumsum(stats2.time),'-o','Color',cm(2,:));
xlabel('iteration');
title('time (s)');
axis tight; padAxis();

% deltaB
subplot(2,2,3);
semilogy(1:nIters,stats1.deltaB,'-o','Color',cm(1,:)); hold on;
semilogy(1:nIters,stats2.deltaB,'-o','Color',cm(2,:));
xlabel('iteration');
title('deltaB');
axis tight; padAxis();

% deltaD
subplot(2,2,4);
semilogy(1:nIters,stats1.deltaD,'-o','Color',cm(1,:)); hold on;
semilogy(1:nIters,stats2.deltaD,'-o','Color',cm(2,:));
xlabel('iteration');
title('deltaD');
axis tight; padAxis();
%--------------------------------------------------------------------------
