%% Fast solver for \|b - Af\|^2 + \lambda \|f\|^2

rng(42);

m = 128;
n = 256;
lambda = 1;

Dm = spdiags([-ones(m,1), ones(m,1)],[0, 1],m,m);
Dm(m,1) = 1;

Dn = spdiags([-ones(n,1), ones(n,1)],[0, 1],n,n);
Dn(n,1) = 1;

A = [kron(Dn, speye(m));
     kron(speye(n), Dm)];

b = randn(2 * m * n,1);

% Least-squares
tic;
AA = [A; sqrt(lambda) * speye(m * n)];
bb = [b; zeros(m * n,1)];
x0 = lsqr(AA,bb,[],1000);
toc

% Matrix
tic;
x1 = (lambda * speye(m * n) + A' * A) \ (A' * b);
toc

% FFT
tic;
Ak = fft2(full(reshape(A' * A(:,1),[m, n])));
Hk = 1 ./ (lambda + Ak);
x2 = reshape(ifft2(Hk .* fft2(reshape(A' * b,[m, n]))),[],1);
toc

% Errors
NRMSE1 = norm(x1 - x0) / norm(x0) %#ok
NRMSE2 = norm(x2 - x0) / norm(x0) %#ok

%%

load('cat.mat');
[m, n, d] = size(I);

I = reshape(I,m * n,d)'; % d x mn
M = repmat(M(:)',d,1); % d x mn
L = L'; % d x 3

% Method 1
N1 = pinv(L) * I;

% Method 2
N2 = (L' * L) \ L' * I;

% Method 3
N3 = zeros(3,m * n);
for j = 1:(m * n)
    if any(M(:,j))
        Mj = diag(M(:,j));
        N3(:,j) = (L' * Mj * L) \ L' * Mj * I(:,j);
    end
end

norm(N1(:) - N2(:))
norm(N1(:) - N3(:))
