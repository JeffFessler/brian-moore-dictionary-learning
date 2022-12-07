%
% dictUpdate demo
%
% Brian Moore
% brimoor@umich.edu
%
% February 1, 2017
%

rng(42);

% Model knobs
d = 20;
m = 10;
n = 100;
SNR = 40;
p = 0.05;

% dictUpdate knobs
nIters = 100;
dr     = 2;
ddim   = [4, 5];
tau    = 0.9; % in [0, 1]
flag   = 0;

% Generate data
Dtrue = randn(d,m);
if dr < min(ddim)
    % Low-rank atoms
    for k = 1:m
        [Uk, Sk, Vk] = svd(reshape(Dtrue(:,k),ddim),'econ');
        Dtruek = Uk(:,1:dr) * Sk(1:dr,1:dr) * Vk(:,1:dr)';
        Dtrue(:,k) = Dtruek(:);
    end
end
Dtrue = bsxfun(@rdivide,Dtrue,sqrt(sum(abs(Dtrue).^2,1)));
B = randn(m,n) .* (rand(m,n) < p);
P = corrupt(Dtrue * B,SNR);

% Standard dictUpdate
opts1 = struct();
opts1.nIters = nIters;
opts1.dr     = dr;
opts1.ddim   = ddim;
opts1.Dtrue  = Dtrue;
opts1.accel  = false;
opts1.tau    = (tau + 1) / norm(B)^2;
opts1.flag   = flag;
[~, stats1] = dictUpdate(P,B,opts1);

% Accelerated dictUpdate
opts2 = struct();
opts2.nIters = nIters;
opts2.dr     = dr;
opts2.ddim   = ddim;
opts2.Dtrue  = Dtrue;
opts2.accel  = true;
opts2.tau    = tau / norm(B)^2;
opts2.flag   = flag;
[~, stats2]  = dictUpdate(P,B,opts2);

% Single step
opts3 = struct();
opts3.nIters = -1;
opts3.dr     = dr;
opts3.ddim   = ddim;
opts3.Dtrue  = Dtrue;
opts3.flag   = flag;
[~, stats3]  = dictUpdate(P,B,opts3);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([819, 299]);
cm = linspecer(2);

% cost
subplot(1,2,1);
phndl = zeros(1,3);
phndl(1) = semilogy(1:nIters,stats1.cost,'-o','Color',cm(1,:)); hold on;
phndl(2) = semilogy(1:nIters,stats2.cost,'-o','Color',cm(2,:));
phndl(3) = semilogy([1, nIters],stats3.cost * [1, 1],'k--');
xlabel('iteration');
title('cost');
legend(phndl,'standard','accel','single step');
axis tight; padAxis();

% nrmse
subplot(1,2,2);
phndl = zeros(1,3);
phndl(1) = plot(1:nIters,stats1.nrmse,'-o','Color',cm(1,:)); hold on;
phndl(2) = plot(1:nIters,stats2.nrmse,'-o','Color',cm(2,:));
phndl(3) = plot([1, nIters],stats3.nrmse * [1, 1],'k--');
xlabel('iteration');
title('nrmse');
legend(phndl,'standard','accel','single step');
axis tight; padAxis();
%--------------------------------------------------------------------------
