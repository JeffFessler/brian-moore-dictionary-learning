function Yhat = SVTp_LLR(Y,p,lambda,nb)
% Syntax:   Yhat = SVTp_LLR(Y,p,lambda,nb);
%
% Locally low-rank (LLR) ell-p SVT
%

% Compute sizes
[m n q] = size(Y);
Nb = prod(nb);
Nx = floor(m / nb(1));
Ny = floor(n / nb(2));

% Loop over blocks
Yhat = zeros(m,n,q);
for i = 1:Nx
    for j = 1:Ny
        % Blockwise ell-p SVT
        ii = ((i - 1) * nb(1) + 1):(i * nb(1));
        jj = ((j - 1) * nb(2) + 1):(j * nb(2));
        Yij = reshape(Y(ii,jj,:),[Nb q]);
        Yhat(ii,jj,:) = reshape(SVTp(Yij,p,lambda),[nb q]);
    end
end
