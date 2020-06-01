%%  Reconstruct simulated images 

% Kevin H. Leung
% 12/16/2018

load('/data/kleung8/system_matrix/system_matrix_inverse_128_sparse.mat'); 
load('/data/kleung8/PET_CNN_Segmentation/training_data/sim_pat_lesion_pairs_v1.mat'); 

lesion_image = double(lesion_image); 
lesion_label = double(lesion_label); 
background_image = double(background_image); 

sim_img_1 = zeros(size(lesion_image,1),128,128); 
sim_img_2 = zeros(size(lesion_image,1),128,128); 
sim_img_3 = zeros(size(lesion_image,1),128,128); 
sim_label = zeros(size(lesion_image,1),128,128); 
for i = 1:size(lesion_image,1)
    fprintf('Simulating image #%d\n',i); 
    % Simulate patient image with tumor
    bg_img = squeeze(background_image(i,:,:)); 
    lesion_img = squeeze(lesion_image(i,:,:)); 
    [recon1, recon2, recon3] = simulate_pet_img(bg_img, lesion_img, H);
    sim_img_1(i,:,:) = recon1;
    sim_img_2(i,:,:) = recon2;
    sim_img_3(i,:,:) = recon3; 
    sim_label(i,:,:) = lesion_label(i,:,:);
end

sim_img_1 = single(sim_img_1); 
sim_img_2 = single(sim_img_2); 
sim_img_3 = single(sim_img_3); 
sim_label = single(sim_label); 

save('sim_pat_16_sub_v1.mat','sim_img_1','sim_img_2','sim_img_3','sim_label'); 

function [recon1, recon2, recon3] = simulate_pet_img(background_img, lesion_img, H)
% Generate simulated image by forward projecting both the image and lesion
% and summing in the projection space. Reconstruct image via OSEM (16
% subsets, 3 iterations) 

% Generate projections for each slice
alpha = 11000; % (1023+1167)/2, scale projections to ~400,000 counts
img = background_img .* alpha; % scale image to get reasonable sinogram count
lesion = lesion_img .* alpha;
projection = H * img(:); % forward projection of image
proj_tumor = H * lesion(:); % forward projection of tumor
projection = projection + proj_tumor; % sum of projections

% Keep the noise level consistent across slices (avg 400000 counts)
while sum(projection(:)) > 500000 || sum(projection(:)) < 300000
    if sum(projection(:)) > 500000
        img = img .* 0.75; % adjust alpha
        lesion = lesion .* 0.75;
        projection = H * img(:);
        proj_tumor = H * lesion(:); % forward projection of tumor
        projection = projection + proj_tumor; % sum of projections
    elseif sum(projection(:)) < 300000
        img = img .* 1.25; % adjust alpha
        lesion = lesion .* 1.25;
        projection = H * img(:);
        proj_tumor = H * lesion(:); % forward projection of tumor
        projection = projection + proj_tumor; % sum of projections
    end
end

projection = poissrnd(projection); % add Poisson noise to sinogrami

% Run 3 iterations of OSEM with 16 subsets
recon_img_iter = ones(128,128); % will serve as initial estimate for first iteration
for iter = 1:3
    F = osem_pet(projection, H, recon_img_iter, 16, 128, 1); % 16 subsets
    recon_img_iter = reshape(F,128,128);
    if iter == 1
        recon1 = recon_img_iter; 
    elseif iter == 2
        recon2 = recon_img_iter; 
    elseif iter == 3
        recon3 = recon_img_iter; 
    end

    % Save each iteration
    fprintf('OSEM Reconstruction, Iteration: %d\n',iter); 
end

end
