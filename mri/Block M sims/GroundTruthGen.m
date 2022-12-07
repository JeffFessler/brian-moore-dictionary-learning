%%
% Generate ground truth Block-M low-rank + sparse (+ noise) image
%

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
% Figure sizes
nx = 128;
ny = 128;
nt = 70; % number of frames

% Background image
Ltrue = imresize(im2double(imread('bwlogo.png')),[nx ny]);

% Moving rectangle params
w = 31; % odd
h = 15; % odd
val = 0.25; % 0 = black, 1 = white

% AWGN params (should be near zero for ground truth)
sigma = 0;
%--------------------------------------------------------------------------

% Get rectangle centers
load('ginput_xy.mat');
%[xc yc] = ginput(nt);
%save('ginput_xy.mat','xc','yc');
xc = round(ny * xc);
yc = round(nx * yc);

% Construct ground truth image
dw = (w - 1) / 2;
dh = (h - 1) / 2;
Mtrue = repmat(Ltrue,[1 1 nt]);
for i = 1:nt
    [rx ry] = meshgrid((-dw:dw) + xc(i),(-dh:dh) + yc(i));
    for j = 1:w
        for k = 1:h
            if ((rx(k,j) >= 1) && (rx(k,j) <= ny) && ...
                (ry(k,j) >= 1) && (ry(k,j) <= nx))
                Mtrue(ry(k,j),rx(k,j),i) = val;
            end
        end
    end
    
    % Add noise
    Mtrue(:,:,i) = Mtrue(:,:,i) + sigma * randn(nx,ny);
end
Mtrue = Mtrue - min(Mtrue(:));
Mtrue = Mtrue / max(Mtrue(:)); % Scale data to [0 1]

% Save data
save('BlockM_groundtruth.mat','Mtrue');
