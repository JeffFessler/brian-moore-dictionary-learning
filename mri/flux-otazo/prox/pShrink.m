function Yhat = pShrink(Y,p,lambda,majorize)
% Syntax:   Yhat = pShrink(Y,p,lambda);
%           Yhat = pShrink(Y,p,lambda,majorize);

% Constants
NUM_ITERS = 200; % # of Newton iterations

% Parse inputs
if ~exist('majorize','var') || isempty(majorize)
    % Default majorization flag
    majorize = false;
end

% Compute shrinkage
absY = abs(Y);
if p == 0
    % Hard thresholding
    Yhat = absY .* (absY >= lambda);
elseif p == 1
    % Soft thresholding
    Yhat = max(0,absY - lambda);
elseif majorize
    % Compute the proximal operator for a reasonable majorizer of
    % the ell-p norm
    Yhat = max(0,absY - lambda * absY.^(p - 1));
else
    % Attempt to exactly compute the ell-p proximal solution via
    % Newton's method
    Jp   = @(x,y) x - y + lambda * p * abs(x).^(p - 1);
    Jpp  = @(x) 1 + (lambda * p * (p - 1)) * abs(x).^(p - 2);
    Yhat = absY;
    for i = 1:NUM_ITERS
        Yhat = Yhat - Jp(Yhat,absY) ./ Jpp(Yhat);
        Yhat((Yhat < 0) | isnan(Yhat)) = 0;
    end
end
Yhat = sign(Y) .* Yhat;
