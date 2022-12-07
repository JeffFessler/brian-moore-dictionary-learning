function  A = Emat_xyt(mask,b1)
%A = Emat_xyt(mask,b1)
%
%
%	Implementation of parallel MRI encoding matrix for dynamic MRI data
%	
%	input:
%			mask : ky-kx-t sampling mask (Ny,Nx,Nt)
%           b1 : coil sensitivity maps (Ny,Nx,Nc)
%
%	output: the operator
%
%	(c) Ricardo Otazo 2008

A.adjoint = false;
A.mask = mask;
A.b1 = permute(b1,[1 2 4 3]);
A = class(A,'Emat_xyt');
