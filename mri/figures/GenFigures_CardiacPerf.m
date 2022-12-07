%% Load the data

load('LS_cardiac_ist.mat');
List = abs(Ltrue);
Sist = abs(Strue);
Mist = abs(Ltrue + Strue);

load('LS_cardiac_SLRN.mat');
Lslrn = abs(Ltrue);
Sslrn = abs(Strue);
Mslrn = abs(Ltrue + Strue);

%% Generate movie from data

% Knobs
ds = 1; % spacing
fps = 5; % frames per second
outpath = 'C:\Users\brimoor\Documents\Google Drive\Winter 2014\ISMRM 2014\movies\svt_optshrink_labels_x2'; % output path
%outpath = '/Users/Brian/Google Drive/Winter 2014/ISMRM 2014/movies/svt_optshrink_labels'; % output path
scale = 2; % resizing scale

% Scale data, if requested
if (scale ~= 1)
    % Scale data
    List = resizeFrames(List,scale);
    Sist = resizeFrames(Sist,scale);
    Mist = resizeFrames(Mist,scale);
    Lslrn = resizeFrames(Lslrn,scale);
    Sslrn = resizeFrames(Sslrn,scale);
    Mslrn = resizeFrames(Mslrn,scale);
end
% Organize data
white = 1e6;
[nx ny nt] = size(Mist);
frames = (1:nt);
X = [];
for frame = frames
    Xf = [Mist(:,:,frame),white * ones(ny,ds),List(:,:,frame),white * ones(ny,ds),Sist(:,:,frame);
          white * ones(ds,3 * nx + 2 * ds);
          Mslrn(:,:,frame),white * ones(ny,ds),Lslrn(:,:,frame),white * ones(ny,ds),Sslrn(:,:,frame)];
    X = cat(3,X,Xf);
end
mval = max(vec(X(X ~= white)));
X(X == white) = mval;
X = double(X);
X = X - min(X(:));
X = X / max(X(:));

% Generate movie
%GenerateMovie(X,outpath,fps);
GenerateMovie_labels(X,outpath,fps);

%% Save the frames separately

% Knobs
ds = 1; % spacing
format = 'tiff'; % image format

[nx ny nt] = size(Mist);
%frames = (1:nt);
frames = [2 14 40];
for frame = frames
    X = [Mist(:,:,frame)   List(:,:,frame)  Sist(:,:,frame);
         Mslrn(:,:,frame) Lslrn(:,:,frame) Sslrn(:,:,frame)];
    white = max(X(:));
    X = [Mist(:,:,frame),white * ones(ny,ds),List(:,:,frame),white * ones(ny,ds),Sist(:,:,frame);
         white * ones(ds,3 * nx + 2 * ds);
         Mslrn(:,:,frame),white * ones(ny,ds),Lslrn(:,:,frame),white * ones(ny,ds),Sslrn(:,:,frame)];
    X = X - min(X(:));
    X = X / max(X(:)); % scale to [0 1]
    X = double(X); % make sure image is double
    %X = im2uint16(X);
    
    %imwrite(X,['./frame' num2str(frame) '.' format],format);
    %imwrite(X,['./frame' num2str(frame) '.jpg'],'jpg','BitDepth',16,'Quality',100,'Mode','lossless');
    %imwrite(X,['./frame' num2str(frame) '.jpg'],'jpg');
    imwrite(X,['./frame' num2str(frame) '.png'],'png');
end

%% Save the frames by algorithm

% Knobs
ds = 2; % spacing

FontSize = 14;

[nx ny nt] = size(Mist);
frames = [2 14 28 40];
Nframes = length(frames);

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
%sp1 = '                                             ';
%sp2 = '                                        ';
sp1 = '                          ';
sp2 = '                     ';

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

%export_fig -pdf -transparent ISTframes_cardiac
%print(gcf,'-depsc2','ISTframes_cardiac.eps');

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot SLRN images
%--------------------------------------------------------------------------
%sp1 = '                                             ';
%sp2 = '                                        ';
sp1 = '                          ';
sp2 = '                     ';

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

%export_fig -pdf -transparent SLRNframes_cardiac
%print(gcf,'-depsc2','SLRNframes_cardiac.eps');
%--------------------------------------------------------------------------
