%% Load the data

load('LS_invivo_ist.mat');
List = abs(Ltrue);
Sist = abs(Strue);
Mist = abs(Ltrue + Strue);

load('LS_invivo_SLRN.mat');
Lslrn = abs(Ltrue);
Sslrn = abs(Strue);
Mslrn = abs(Ltrue + Strue);

%% Save the frames by algorithm

% Knobs
ds = 2; % spacing
ROI = [45 160 1 90];
FontSize = 14;

frames = [6 13 20 50];
Nframes = length(frames);

% Restrict images to ROI
List = List(ROI(1):ROI(2),ROI(3):ROI(4),:);
Sist = Sist(ROI(1):ROI(2),ROI(3):ROI(4),:);
Mist = Mist(ROI(1):ROI(2),ROI(3):ROI(4),:);
Lslrn = Lslrn(ROI(1):ROI(2),ROI(3):ROI(4),:);
Sslrn = Sslrn(ROI(1):ROI(2),ROI(3):ROI(4),:);
Mslrn = Mslrn(ROI(1):ROI(2),ROI(3):ROI(4),:);

[ny nx nt] = size(Mist);

w1  = max([max(vec(Mist(:,:,frames)))  max(vec(List(:,:,frames)))  max(vec(Sist(:,:,frames)))]);
w2 = max([max(vec(Mslrn(:,:,frames))) max(vec(Lslrn(:,:,frames))) max(vec(Sslrn(:,:,frames)))]);

Xist = [];
Xslrn = [];
for i = 1:Nframes
    frame = frames(i);
    
    Xist  = [Xist; [Mist(:,:,frame) ,w1 * ones(ny,ds),List(:,:,frame), w1 * ones(ny,ds),Sist(:,:,frame)]]; %#ok
    Xslrn = [Xslrn;[Mslrn(:,:,frame),w2 * ones(ny,ds),Lslrn(:,:,frame),w2 * ones(ny,ds),Sslrn(:,:,frame)]]; %#ok
    
    if (i < length(frames))
        Xist  = [Xist; w1 * ones(ds,3 * nx + 2 * ds)]; %#ok
        Xslrn = [Xslrn;w2 * ones(ds,3 * nx + 2 * ds)]; %#ok
    end
end

Xist = Xist - min(Xist(:));
Xist = Xist / max(Xist(:)); % scale to [0 1]
Xist = double(Xist); % make sure image is double
%Xist = im2uint16(Xist);

Xslrn = Xslrn - min(Xslrn(:));
Xslrn = Xslrn / max(Xslrn(:)); % scale to [0 1]
Xslrn = double(Xslrn); % make sure image is double
%Xslrn = im2uint16(Xslrn);

%--------------------------------------------------------------------------
% Plot IST images
%--------------------------------------------------------------------------
sp1 = '                  ';
sp2 = '                   ';

figure;
imshow(Xist);

title(['L + S' sp1(1:(end-5)) 'L' sp1 'S  '],'FontSize',FontSize);

ystr = '';
for i = Nframes:-1:2
    ystr = [ystr 'F' num2str(frames(i)) sp2]; %#ok
end
ystr = [ystr 'F' num2str(frames(1))];
yh = ylabel(ystr,'FontSize',FontSize);
set(yh, 'Units', 'Normalized');
pos = get(yh, 'Position');
set(yh, 'Position', pos + [0.085 0 0]);

TightFigure(gca);

%print(gcf,'-depsc2','ISTframes_invivo.eps');

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot SLRN images
%--------------------------------------------------------------------------
sp1 = '                  ';
sp2 = '                   ';

figure;
imshow(Xslrn);

title(['L + S' sp1(1:(end-5)) 'L' sp1 'S  '],'FontSize',FontSize);

ystr = '';
for i = Nframes:-1:2
    ystr = [ystr 'F' num2str(frames(i)) sp2]; %#ok
end
ystr = [ystr 'F' num2str(frames(1))];
yh = ylabel(ystr,'FontSize',FontSize);
set(yh, 'Units', 'Normalized');
pos = get(yh, 'Position');
set(yh, 'Position', pos + [0.085 0 0]);

TightFigure(gca);

%print(gcf,'-depsc2','SLRNframes_invivo.eps');
%--------------------------------------------------------------------------
