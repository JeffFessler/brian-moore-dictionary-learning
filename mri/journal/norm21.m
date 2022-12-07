function nY = norm21(Y)
% Syntax:   nY = norm21(Y);
%
% Mixed 2-1 norm
%

% Columns
%nY = sum(sqrt(sum(Y.^2,1)));

% Rows
nY = sum(sqrt(sum(Y.^2,2)));
