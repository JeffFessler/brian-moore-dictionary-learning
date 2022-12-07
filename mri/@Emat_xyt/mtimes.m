function y = mtimes(A,x)

% Perform multiplication
if (A.adjoint == true)
    % y = A' * x
    arg0 = ifft2c_mri(bsxfun(@times,x,A.mask));
    arg1 = sum(bsxfun(@times,arg0,conj(A.b1)),4);
    arg2 = squeeze(sum(abs(A.b1).^2,4));
    y = bsxfun(@rdivide,arg1,arg2);
else
    % y = A * x
    y = bsxfun(@times,fft2c_mri(bsxfun(@times,x,A.b1)),A.mask);
end
