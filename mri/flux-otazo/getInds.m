function inds = getInds(idx,N,MAX_IDX)
% Syntax:   inds = getInds(idx,N,MAX_IDX);

% Get indices
Navg = N / MAX_IDX;
Nceil = ceil(Navg);
Nfloor = floor(Navg);
N1 = ceil(MAX_IDX * (Navg - Nfloor));
if (idx <= N1)
    inds = ((idx - 1) * Nceil + 1):min(N,idx * Nceil);
else
    Ninc = N1 * Nceil;
    idx2 = idx - N1;
    inds = Ninc + (((idx2 - 1) * Nfloor + 1):min(idx2 * Nfloor,N - Ninc));
end

% Permute indices around search cube
rng(0);
Nperm = randperm(N);
inds = Nperm(inds);
