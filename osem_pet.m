function [ F ] = osem_pet(projection, H, initial_estimate, subsets, angle_bin, num_iter)
%OSEM Algorithm
%   Implementation of OSEM for PET reconstruction
%   Kevin H. Leung
%   11/30/17
% 
%   F - Reconstructed image (dim - n) (save each iteration) 
%   H - System matrix (dim - m,n)
%   G - Projection (measurement) (dim - m)
%   projection - should be given to the function as an unrolled vector
%   s - point sensitivity vector (dim - n)
%   initial_estimate - initial estimate of image (usually matrix of ones)
%   subsets - number of subsets 
%   angle_bin - number of angular bins 
%   num_iter - number of full iterations (num of subiters = num_subsets*num_iter)

G = projection; 
G = reshape(G, length(G)/angle_bin, angle_bin);
F = initial_estimate(:);
for i = 1:num_iter
    for j = 1:subsets
        subset = j:subsets:angle_bin; % get different subset of projection data
        G_sub = G(:,subset);
        tmp = reshape(H*F, length(G(:))/angle_bin, angle_bin);
        tmp_sub = tmp(:,subset);
        for n = 1:length(F) % update every element in F (reconstruction)
            H_sub = reshape(H(:,n), length(G(:))/angle_bin, angle_bin); 
            H_sub = H_sub(:, subset); 
            s = sum(H_sub(:)); % nth element of sensitivity vector for current subset 
            F(n) =  F(n) * (1/s) * sum((G_sub(:) ./ (tmp_sub(:)+eps)) .* H_sub(:));
        end
    end
end

end

