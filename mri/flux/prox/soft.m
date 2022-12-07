function Yhat = soft(Y,lambda)
% Syntax:   Yhat = soft(Y,lambda);
%
% Soft-thresholding
%

% Apply soft-thresholding
Yhat = max(abs(Y) - lambda,0) .* sign(Y);
Yhat(isnan(Yhat)) = 0;
