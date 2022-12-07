function [Ys mask] = subsampleData(Y,p,ps)
% Syntax:   [Ys mask] = subsampleData(Y,p,ps);

% Parse inputs
[ny nx nt nc] = size(Y);
num = round(p * nx);

% Undersample data
Ys = zeros(ny,nx,nt,nc);
mask = false(ny,nx,nt);
for i = 1:nt
    cols = sampleCols(num,ps);
    mask(:,cols,i) = true;
    Ys(:,cols,i,:) = Y(:,cols,i,:);
end
