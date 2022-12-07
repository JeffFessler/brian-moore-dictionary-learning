function Yhat = FroShrink(Y,lambda)
% Syntax:   Yhat = FroShrink(Y,lambda);
%
% Frobenius shrinkage operator
%

% Perform Frobenius shrinkage
Yhat = (1 / (1 + lambda)) * Y;
