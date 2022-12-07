function Yhat = SVT_LLR(Y,lambda,nb)
% Syntax:   Yhat = SVT_LLR(Y,lambda,nb);
%
% Locally low-rank (LLR) singular value thresholding
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
        % Perform singular value thresholding
        ii = ((i - 1) * nb(1) + 1):(i * nb(1));
        jj = ((j - 1) * nb(2) + 1):(j * nb(2));
        Yij = reshape(Y(ii,jj,:),[Nb q]);
        Yhat(ii,jj,:) = reshape(SVT(Yij,lambda),[nb q]);
    end
end
