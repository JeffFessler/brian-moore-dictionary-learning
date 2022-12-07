function inds = sampleCols(num,ps)
% Syntax:   inds = sampleCols(num,ps);

% Sample columns based on empirical sampling frequencies, ps
inds = [];
while numel(inds) < num
    P       = cumsum(ps);
    [~,idx] = histc(rand(1,num),[0, P / P(end)]);
    idx     = unique(idx);
    ps(idx) = 0;
    inds    = [inds, idx(randperm(numel(idx)))]; %#ok
end
inds = sort(inds(1:num));
