function FX = ifft2c_mri(X)

FX = fftshift(fft(fftshift(X,1),[],1),1) / sqrt(size(X,1));
FX = fftshift(fft(fftshift(FX,2),[],2),2) / sqrt(size(X,2));
