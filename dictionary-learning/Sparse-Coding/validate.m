%
% sparseCoding validation
%
% Brian Moore
% brimoor@umich.edu
%
% January 31, 2017
%

rng(42);

% Model knobs
d = 10;
m = 10;
n = 20;

% sparseCoding knobs
mu     = 1;
nIters = 50;
accel  = true;
tau    = 0.99;
flag   = 0;

% Generate problem instance
P = randn(d,n);
D = randn(d,m);

% Requires CVX from here: http://cvxr.com/cvx/download/
tic;
cost = @(B) 0.5 * square_pos(norm(P - D * B,'fro')) + mu * norm(vec(B),1);
cvx_begin quiet
    variables Bcvx(m,n)
    minimize(cost(Bcvx));
cvx_end
cost_cvx = cost(Bcvx);
time_cvx = toc %#ok

% sparseCoding
tic;
opts.nIters  = nIters;
opts.type    = 'soft';
opts.accel   = accel;
opts.tau     = tau / norm(D)^2;
opts.flag    = flag;
[~, stats]   = sparseCoding(P,D,mu,opts);
time_pg      = toc %#ok

% Plot results
figure;
phndl = zeros(1,2);
phndl(1) = semilogy(1:nIters,stats.cost, 'b-o'); hold on;
phndl(2) = semilogy([1, nIters],cost_cvx * [1, 1],'k--');
xlabel('iteration');
ylabel('cost');
legend(phndl,'sparseCoding','CVX');
axis tight; padAxis();
