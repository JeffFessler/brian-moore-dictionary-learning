function full_data = data_share_fill(undersamp_data,sampling_pattern)
%function full_data = data_share_fill(undersamp_data,sampling_pattern)
% data sharing
% takes undersampled dynamic data and uses data-sharing (0th order interp) to fill in gaps
% output is "fully sampled"
% data-sharing fills in unsampled k-space locations with values from nearest (in time) sampling of that location, using mean for ties
% leaves zeros if k-space location not sampled
%
% inputs:
%	sampling_pattern [Nx Ny Nf] logical array
%	undersamp_data [Ns Nc] where Ns = a*Nx*Ny*Nf and a is undersamp factor
%
% outputs:
%	full_data [Nx Ny Nf Nc] complex array
% 02/05/14 Mai Le


Nx = size(sampling_pattern,1);
Ny = size(sampling_pattern,2);
Nf = size(sampling_pattern,3);

Ns = size(undersamp_data,1);
Nc = size(undersamp_data,2); 

flip_samp = flipdim(sampling_pattern,3);

samps = {sampling_pattern,flip_samp};

zfill = reshape(embed(undersamp_data, logical(sampling_pattern)),Nx,Ny,Nf,Nc);

for ii = 1:2
	curr_samp = samps{ii};
	%interp_ind = cumsum(int32(curr_samp),3);
	interp_ind = cumsum((curr_samp),3);
	%interp_ind_pos = interp_ind + int32(interp_ind == 0);
	interp_ind_pos = interp_ind + (interp_ind == 0);
	siip(:,:,:,ii) = interp_ind_pos; % keep!!
    
	%full_interp_ind_pos = Nx*Ny*(interp_ind_pos(:)-1) + int32(repmat((0:Nx*Ny-1)', [Nf 1])) + 1;
	full_interp_ind_pos = Nx*Ny*(interp_ind_pos(:)-1) + (repmat((0:Nx*Ny-1)', [Nf 1])) + 1;
    
	if Nc > 1
		%full_interp_ind_pos = repmat(full_interp_ind_pos,[Nc 1]) + int32(kron((0:Nc-1)',Nx*Ny*Nf*ones(Nx*Ny*Nf,1)));	
		full_interp_ind_pos = repmat(full_interp_ind_pos,[Nc 1]) + (kron((0:Nc-1)',Nx*Ny*Nf*ones(Nx*Ny*Nf,1)));	
    end
    
	[trash, inds] = sort(curr_samp,3,'descend');
	sinds(:,:,:,ii) = inds; % keep!!
    
	full_inds = Nx*Ny*(inds(:)-1) + repmat((1:Nx*Ny)',[Nf 1]);
	if Nc > 1
		full_inds = repmat(full_inds,[Nc 1]) + kron((0:Nc-1)',Nx*Ny*Nf*ones(Nx*Ny*Nf,1));
    end
    
	if ii == 2
		curr_zfill = vec(flipdim(zfill,3));
	else
		curr_zfill = vec(zfill);
	end
	nonzero_elms(:,ii) = curr_zfill(full_inds);
	
	zero_less_fill = reshape(nonzero_elms(:,ii),Nx,Ny,Nf,Nc);


	curr_full_data = reshape(zero_less_fill(full_interp_ind_pos),Nx,Ny,Nf,Nc);

	if ii == 2 
		curr_full_data  = flipdim(curr_full_data ,3);
	end

	sfull_data(:,:,:,:,ii) = curr_full_data;
end
siip1 = double(siip(:,:,:,1));
c1_ndx = (siip1-1)*Nx*Ny+repmat(reshape(0:Nx*Ny-1,Nx,Ny),[1 1 Nf])+1;
frames1 = vec(sinds(:,:,:,1));
c1 = vec(reshape(frames1(c1_ndx),Nx,Ny,Nf));


siip2 = double(flipdim(siip(:,:,:,2),3));
siip2 = repmat(siip2(:,:,1) + 1,[1 1 Nf]) - siip2;
c2_ndx = (siip2-1)*Nx*Ny+repmat(reshape(0:Nx*Ny-1,Nx,Ny),[1 1 Nf])+1;
frames2 = vec(sinds(:,:,:,1)); % need to have 1!
c2 = vec(reshape(frames2(c2_ndx),Nx,Ny,Nf));



c3 = vec(reshape(kron((1:Nf),ones(Nx*Ny,1)),Nx,Ny,Nf));
%[c1 c2 c3 vec(1:Nx*Ny*Nf)]
choose_flip = repmat(reshape(abs(c2 - c3) < abs(c1 - c3),Nx,Ny,Nf),[1 1 1 Nc]);

tie = repmat(reshape((abs(c2 - c3) == abs(c1 - c3)) & (c1 ~= c2)  ,Nx,Ny,Nf),[1 1 1 Nc]); 

flip_weights = choose_flip + 0.5*tie;
full_data = sfull_data(:,:,:,:,1).*(1-flip_weights) + sfull_data(:,:,:,:,2).*flip_weights;

% Vectorize data
function vx = vec(x)
vx = x(:);
