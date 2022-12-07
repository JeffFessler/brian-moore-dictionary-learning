function Y = resizeFrames(X,arg)
% Syntax:   Y = resizeFrames(X,scale);
%           Y = resizeFrames(X,[nrows ncols]);

% Resize frames
[ny nx] = size(imresize(X(:,:,1),arg));
nt = size(X,3);
Y = zeros(ny,nx,nt);
for i = 1:nt
    Y(:,:,i) = imresize(X(:,:,i),arg);
end
