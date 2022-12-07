function Yhat = OptShrink_LLR(Y,r,nb)
% Syntax:   Yhat = OptShrink_LLR(Y,r,nb);
%
% Locally low-rank (LLR) OptShrink
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
        % Blockwise OptShrink
        ii = ((i - 1) * nb(1) + 1):(i * nb(1));
        jj = ((j - 1) * nb(2) + 1):(j * nb(2));
        Yij = reshape(Y(ii,jj,:),[Nb q]);
        Yhat(ii,jj,:) = reshape(OptShrink(Yij,r),[nb q]);
    end
end
