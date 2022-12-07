%% Generate SENSE/SAMP images

erase;
load('cardiac_perf_R8.mat');

for idx = 1:8
    %I = abs(b1(:,:,idx));
    I = double(abs(kdata(:,:,idx,1))' > 0);
    
    I = I - min(I(:));
    I = I / max(I(:));
    I = round(255 * I);
    
    scrsz = get(0,'Screensize');
    xyc = 0.5 * scrsz(3:4);
    dim = [size(b1,2) size(b1,1)];
    
    f = figure('Menubar','none', ...
               'Position',[(xyc - 0.5 * dim) dim]);
    ax = axes('Position',[0 0 1 1]);
    image([1 dim(1)],[1 dim(2)],I,'Parent',ax);
    colormap(gray(256));
    axis(ax,'off');
    drawnow; pause(0.01);
    
    %path = sprintf('/Users/Brian/Desktop/MRI/sense%i',idx);
    path = sprintf('/Users/Brian/Desktop/MRI/samp%i',idx);
    export_fig('-pdf','-transparent',path);
    drawnow; pause(0.01);
    delete(f);
end
