function [blocks,idx, rows, cols, temp] = my_im2col_3D(I,blkSize, slidingDis, tSlideingDis)
[aa, bb, cc] = size(I);
idxMat = zeros([aa, bb, cc]- blkSize + 1);
idxMat([[1:slidingDis:end-1],end],[[1:slidingDis:end-1], end],[[1:tSlideingDis:end-1], end]) = 1; % take blocks in distances of 'slidingDix', but always take the first and last one (in each row and column).
idx = find(idxMat);
[rows,cols, temp] = ind2sub(size(idxMat),idx);
blocks = zeros(prod(blkSize),length(idx));
for i = 1:length(idx)
    currBlock = I(rows(i):rows(i)+blkSize(1)-1,cols(i):cols(i)+blkSize(2)-1, temp(i):temp(i)+blkSize(3)-1);
    blocks(:,i) = currBlock(:);
end
