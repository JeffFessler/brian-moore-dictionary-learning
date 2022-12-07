%% Visualize dictionary atoms

% Knobs
inpath = 'otazo_R8_lads_p1_mse250.mat';
pdim   = [8, 8, 5];

% Load data
load(inpath,'Dhat');

% Construct atom images
D  = reshape(Dhat,8,8,5,320);
D  = squeeze(D(:,:,1,:));
Dc = cell(8,8,5);
for i = 1:320
    Dc{i} = abs(D(:,:,i));
    %Dc{i} = real(D(:,:,i));
end

% Display atoms
figure;
for i = 1:5
    subplot(2,3,i);
    Di = cell2mov(Dc(:,:,i),1,nan);
    Di = Di / nanmax(Di(:));
    Di(isnan(Di)) = 1;
    imshow(Di);
end
