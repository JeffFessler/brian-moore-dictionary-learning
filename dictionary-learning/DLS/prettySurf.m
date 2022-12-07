function F = prettySurf(f,M,zlim)
% Syntax:   F = prettySurf(f,M);
%           F = prettySurf(f,M,zlim);

if exist('zlim','var') && ~isempty(zlim)
    f = max(zlim(1),min(f,zlim(2)));
end
F = reshape(f,size(M));
F = F - min(F(M));
F(~M) = nan;
