function X = fft2c_mri(FX)

X = sqrt(size(FX,1)) * fftshift(ifft(fftshift(FX,1),[],1),1);
X = sqrt(size(FX,2)) * fftshift(ifft(fftshift(X,2),[],2),2);
