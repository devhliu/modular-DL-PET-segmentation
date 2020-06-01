%% Extract information about tumor size, shape, and intensity information from patient data 

% Kevin H. Leung
% 12/16/2018

clear; clc; 
close all; 

%% Find relevant patient slices and corresponding manual segmentations

% load existing patients and find the corresponding indicies 
load('C:\Users\khleu\OneDrive\Documents\MATLAB\Research_PET_Segmentation\Segm_methods\scripts\patient_num_info.mat');
pat_num = string(pat_num);

lungData = readtable('C:\Users\khleu\OneDrive\Documents\MATLAB\PET Images\Lung_data_sparse.xlsx'); 
lungData = table2cell(lungData); 
lungData = lungData(:,1:5); 
pat_no = char(lungData(:,1));

pat_no_str = cell(length(pat_no),1); 
for i = 1:length(pat_no)
    if length(strtrim(pat_no(i,:))) < 7 
        pat_no_str{i} = ['0', strtrim(pat_no(i,:))]; 
    elseif length(strtrim(pat_no(i,:))) > 7
        tmp = strtrim(pat_no(i,:)); 
        pat_no_str{i} = tmp(3:end); 
        if length(tmp(3:end)) < 7 
            pat_no_str{i} = ['0', tmp(3:end)]; 
        end
    else
        pat_no_str{i} = strtrim(pat_no(i,:)); 
    end
end

pat_no_str = string(pat_no_str); 

pat_mapping = cell(size(pat_num)); 
for i = 1:length(pat_num)
    for j = 1:length(pat_no_str)
        if pat_num(i)==pat_no_str(j)
            pat_mapping{i} = j; 
            break; 
        end
    end
end

pat_mapping = cell2mat(pat_mapping);

% Re-order mapping of patients 
lungData = lungData(pat_mapping,:); 

 % get subset of patients with Discovery LS scanner 
scanner_subset = (string(lungData(:,3))=="Discovery LS");
scanner_subset_other = (string(lungData(:,3))~="Discovery LS");

other_lungData = lungData(scanner_subset_other,:); 
lungData = lungData(scanner_subset,:); 

% Get new patient numbers with scanner Discovery LS scanner 
pat_no = char(lungData(:,1));
pat_no_str = cell(length(pat_no),1); 
for i = 1:length(pat_no)
    if length(strtrim(pat_no(i,:))) < 7 
        pat_no_str{i} = ['0', strtrim(pat_no(i,:))]; 
    elseif length(strtrim(pat_no(i,:))) > 7
        tmp = strtrim(pat_no(i,:)); 
        pat_no_str{i} = tmp(3:end); 
        if length(tmp(3:end)) < 7 
            pat_no_str{i} = ['0', tmp(3:end)]; 
        end
    else
        pat_no_str{i} = strtrim(pat_no(i,:)); 
    end
end

pat_no_str = char(pat_no_str); 

% Load all patient slices and manual segmentation slices 
PETimg_files = dir('H:\Backup_research_3_26_2018\PET_NGS\Lung_cancer\image_data\pre_therapy\Complete_data\PETimg\*.mat'); 
manSeg_files = dir('H:\Backup_research_3_26_2018\PET_NGS\Lung_cancer\image_data\pre_therapy\Complete_data\Manual_segm\*.mat'); 

PETimg = cell(length(pat_no_str),1); 
manSeg = cell(length(pat_no_str),1); 
tumor_slices = cell(length(pat_no_str),1); % logical index of tumor slices
for i = 1:length(pat_no_str)
    % Get corresponding PET images 
    for j = 1:length(PETimg_files)
        tmp = PETimg_files(j).name;
        tmp = tmp(1:7);
        if strcmp(tmp,pat_no_str(i,:))
            load([PETimg_files(j).folder, '\', PETimg_files(j).name]);
            PETimg{i} = vol_vals;
            break; 
        end
    end
    % Get corresponding manual segmentation 
    for j = 1:length(manSeg_files)
        tmp = manSeg_files(j).name;
        tmp = tmp(1:7);
        if strcmp(tmp,pat_no_str(i,:))
            load([manSeg_files(j).folder, '\', manSeg_files(j).name]);
            manSeg{i} = manual_seg;
            tmp_slice = squeeze(sum(sum(manual_seg))); % find tumor slice indices
            tumor_slices{i} = logical(tmp_slice>0); 
            break; 
        end
    end
end

%% Extract patient dependent information on tumor characteristics 

% Tumor volume (height and width lengths)
tumor_volume = cell(length(tumor_slices),1);
tumor_dimension = cell(length(tumor_slices),1); 
tumor_diameter = cell(length(tumor_slices),1); 

for i = 1:length(tumor_slices)
    idx = tumor_slices{i}; 
    seg = manSeg{i}; 
    seg = seg(:,:,idx); 
    tmp_vol = zeros(size(seg,3),1); 
    tmp_dim = zeros(size(seg,3),2); 
    tmp_dia = zeros(size(seg,3),1); 
    for j = 1:size(seg,3)
        tmp_seg = seg(:,:,j); 
        tmp_vol(j) = sum(tmp_seg(:)); 
        x = find(tmp_seg>0); 
        [I,J] = ind2sub(size(tmp_seg),x); 
        tmp_dim(j,1) = max(I)-min(I); % find row dimension of tumor
        tmp_dim(j,2) = max(J)-min(J); % find col dimension of tumor
        tmp_dia(j) = sqrt(sum(tmp_dim(j,:).^2)); 
    end
    tumor_volume{i} = tmp_vol; 
    tumor_dimension{i}  = tmp_dim; 
    tumor_diameter{i} = tmp_dia; 
end

tumor_vol = []; 
tumor_dim = []; 
tumor_dia = []; 
for i = 1:length(tumor_volume)
    tmp = tumor_volume{i}; 
    tumor_vol = [tumor_vol; tmp]; 
    tmp = tumor_dimension{i}; 
    tumor_dim = [tumor_dim; tmp]; 
    tmp = tumor_diameter{i}; 
    tumor_dia = [tumor_dia; tmp]; 
end

% Tumor intensity and variance per slice 
tumor_intensity = cell(length(tumor_slices),1);
tumor_variance = cell(length(tumor_slices),1); 
for i = 1:length(tumor_slices)
    idx = tumor_slices{i}; 
    seg = manSeg{i};
    img = PETimg{i};
    seg = seg(:,:,idx);
    img = img(:,:,idx); 
    if i == 29 || i == 37
        img = img ./ 10^4; % balance out scaling factor
    end
    tmp_int = zeros(size(seg,3),1); 
    tmp_var = zeros(size(seg,3),1); 
    for j = 1:size(seg,3)
        tmp_seg = seg(:,:,j); 
        tmp_img = img(:,:,j); 
        tmp_tum = tmp_seg .* tmp_img; 
        tmp_nonzero = (tmp_tum(:) > 0); 
        tmp_tum = tmp_tum(:); 
        tmp_int(j) = mean(tmp_tum(tmp_nonzero)); % mean intensity per slice 
        tmp_var(j) = var(tmp_tum(tmp_nonzero));         
    end
    tumor_intensity{i} = tmp_int; 
    tumor_variance{i} = tmp_var; 
end

tumor_int = []; 
tumor_var = []; 
for i = 1:length(tumor_intensity)
    tmp = tumor_intensity{i}; 
    tumor_int = [tumor_int; tmp]; 
    tmp = tumor_variance{i}; 
    tumor_var = [tumor_var; tmp]; 
end

% Patient backgorund intensity around tumor
background_intensity = cell(length(tumor_slices),1);
background_variance = cell(length(tumor_slices),1); 
for i = 1:length(tumor_slices)
    idx = tumor_slices{i}; 
    seg = manSeg{i};
    img = PETimg{i};
    seg = seg(:,:,idx);
    img = img(:,:,idx); 
    if i == 29 || i == 37
        img = img ./ 10^4; % balance out scaling factor
    end
    tmp_int = zeros(size(seg,3),1); 
    tmp_var = zeros(size(seg,3),1); 
    for j = 1:size(seg,3)
        tmp_seg = seg(:,:,j); 
        tmp_img = img(:,:,j); 
        x = find(tmp_seg>0); 
        [row,col] = ind2sub(size(tmp_seg),x); 
        % Circular ROI
        radius = sqrt(length(min(row):max(row))^2 + length(min(col):max(col))^2)/2; % 200% radius of tumor (diameter of true square ROI)
        cir_y = round((min(row)+max(row))/2);
        cir_x = round((min(col)+max(col))/2);
        [xx,yy] = ndgrid((1:128)-cir_y,(1:128)-cir_x);
        cir_img_roi = (xx.^2 + yy.^2)<radius^2; % circular ROI around tumor (200% tumor radius)
        cir_img_roi = cir_img_roi .* tmp_img;
        cir_img_roi = cir_img_roi .* ~tmp_seg; 
        tmp_nonzero = (cir_img_roi(:) > 0); 
        tmp_roi = cir_img_roi(:); 
        tmp_int(j) = mean(tmp_roi(tmp_nonzero)); % mean intensity per slice 
        tmp_var(j) = var(tmp_roi(tmp_nonzero)); 
    end
    background_intensity{i} = tmp_int; 
    background_variance{i} = tmp_var; 
end

background_int = []; 
background_var = []; 
for i = 1:length(background_intensity)
    tmp = background_intensity{i}; 
    background_int = [background_int; tmp]; 
    tmp = background_variance{i}; 
    background_var = [background_var; tmp]; 
end

% Tumor to background ratio 
tumor_background_intensity = cell(length(tumor_intensity),1); 
tumor_bg_int = []; 
for i = 1:length(tumor_intensity) 
    tmp_tumor = tumor_intensity{i}; 
    tmp_background = background_intensity{i}; 
    tmp_tumor_bg = tmp_tumor ./ tmp_background; 
    tumor_background_intensity{i} = tmp_tumor_bg; 
    if tmp_tumor_bg > 1.0
        tumor_bg_int = [tumor_bg_int; tmp_tumor_bg]; 
    end
end

% Tumor shape (elliptical fourier coefficients) using 10 harmonic
% coefficients
eFourierCoef = cell(length(tumor_slices),1);
ct = 0; 
pat_tumor_shape = cell(1); 
true_shape = cell(1); 
for i = 1:length(tumor_slices)
    idx = tumor_slices{i};
    seg = manSeg{i};
    seg = seg(:,:,idx);
    tmp_eFCoef = []; 
    for j = 1:size(seg,3)
        tmp_seg = seg(:,:,j);
        [x, y] = ind2sub(size(tmp_seg), find(tmp_seg>0)); 
        row = length(min(x):max(x)); 
        col = length(min(y):max(y)); 
        diff_row = 41 - (max(x)-min(x)); 
        diff_col = 31 - (max(y)-min(y)); 
        tmp = zeros(41, 31);
        tmp(round(diff_row/2):round(diff_row/2)+row-1, round(diff_col/2):round(diff_col/2)+col-1) = tmp_seg(min(x):max(x), min(y):max(y)); 
        outline = bwboundaries(tmp);
        outline = outline{1};
        rFSDs = fEfourier(outline, 5, 1, 1); % normalize for size and orientation
        outln = rEfourier(rFSDs, 5,100);
        outln = abs(flipud(outln.*50)); % increase resolution
        bw = poly2mask(outln(:,2),outln(:,1),41.*50,31.*50);
        [row, col] = ind2sub(size(bw), find(bw>0));
        
        % Test for valid polygon (no self-intersections)
        [x0,y0,segments]=selfintersect(outln(:,2),outln(:,1));
        if isempty(x0)
            ct = ct + 1;
            tmp_eFCoef = cat(3, tmp_eFCoef, rFSDs);
            pat_tumor_shape{ct} = bw(min(row):max(row), min(col):max(col));
            true_shape{ct} = tmp_seg(min(x):max(x), min(y):max(y)); 
        end
    end
    eFourierCoef{i} = tmp_eFCoef;
end

% Get elliptical fourier coefficients from all patient slices with valid shapes 
eFC_all = zeros(4,5,ct);
ct = 0;
for i = 1:length(eFourierCoef)
    tmp = eFourierCoef{i};
    mini_ct = 0;
    for j = 1:size(tmp,3)
        ct = ct + 1;
        mini_ct = mini_ct + 1;
        eFC_all(:,:,ct) = tmp(:,:,mini_ct);
    end
end

% Tumor statistics are in the following variables: 
% tumor_vol, tumor_dim, tumor_int, tumor_var, background_int,
% background_var, tumor_bg_intensity, eFC_all


%% Simulate tumors in background patient slices 
load('filt_pat_bg_slices_v2.mat');
load('tumor_centers_pat_bg.mat'); 

lesion_image = zeros(1740,128,128);
lesion_label = zeros(1740,128,128); 
background_image = zeros(1740,128,128); 

% Select points for background slices 
ct = 0; 
for i = 1:length(filtered_slices_v2)
    
    idx = filtered_slices_v2{i}; 
    img = PETimg{i}; 
    img = img(:,:,idx); 
    if i == 29 || i == 37
        img = img ./ 10^4; % balance out scaling factor
    end
    
    points = tumor_centers{i}; 
    
    for j = 1:size(img,3)
        ct = ct + 1; 
        bg_img = img(:,:,j); 
        point = points(j,:); 
        
        % Offset the points in background slice by +- 2 
        x = point(1) + round(unifrnd(-2,2)); 
        y = point(2) + round(unifrnd(-2,2)); 
        
        % Generate random lesion 
        if rand < 0.5
            hetero = 0;
        else
            hetero = 1;
        end
        [lesion, mask] = generate_tumor(eFC_all, tumor_dia, tumor_vol, hetero, 100);
        
        % Place lesion 
        height = size(lesion, 1);
        width = size(lesion,2);
        x_ = round(x - width/2);
        y_ = round(y - height/2);
        lesion_img = zeros(128);
        lesion_img(y_:height+y_-1, x_:width+x_-1) = lesion;
        lesion_mask = zeros(128);
        lesion_mask(y_:height+y_-1, x_:width+x_-1) = mask;
        
        % Adjust for tumor-to-background ratio using circular ROI
        [row,col] = ind2sub(size(lesion_mask),find(lesion_mask>0));
        radius = sqrt(length(min(row):max(row))^2 + length(min(col):max(col))^2)/2; % radius of tumor (diameter of true square ROI)
        cir_y = round((min(row)+max(row))/2);
        cir_x = round((min(col)+max(col))/2);
        [xx,yy] = ndgrid((1:128)-cir_y,(1:128)-cir_x);
        cir_img_roi = (xx.^2 + yy.^2)<radius^2; % circular ROI around tumor (200% tumor radius)
        cir_img_roi = cir_img_roi .* bg_img; % get ROI
        mean_bg_intensity = mean(cir_img_roi(cir_img_roi>0)); % mean background intensity
        
        % Model tumor-to-background ratio via Beta distribution
        tumor_to_bg = generate_TBR(tumor_bg_int);
        mean_tumor_intensity = tumor_to_bg * mean_bg_intensity; % - mean_bg_intensity;
        
        lesion_img = generate_tumor_intensity(lesion_img, mean_tumor_intensity, tumor_var, mean_bg_intensity, background_var);
        lesion_image(ct,:,:) = lesion_img;
        lesion_label(ct,:,:) = lesion_mask;
        background_image(ct,:,:) = bg_img;
    end
    
end

lesion_image = single(lesion_image); 
lesion_label = single(lesion_label); 
background_image = single(background_image); 

% Check 
for i = 1:10:size(lesion_image,1)
    subplot(1,3,1); imagesc(squeeze(lesion_image(i,:,:))); 
    subplot(1,3,2); imagesc(squeeze(lesion_label(i,:,:))); 
    subplot(1,3,3); imagesc(squeeze(background_image(i,:,:)) + squeeze(lesion_label(i,:,:))); 
    title(sprintf('%d',i)); 
    if i == 1
        pause; 
    else
        pause(0.25);
    end
end

subplot(1,3,1); imagesc(squeeze(lesion_image(end,:,:)));
subplot(1,3,2); imagesc(squeeze(lesion_label(end,:,:)));
subplot(1,3,3); imagesc(squeeze(background_image(end,:,:)) + squeeze(lesion_label(end,:,:)));

save('sim_pat_lesion_pairs_v7.mat','lesion_image','lesion_label','background_image'); 

%% Helper functions for simulating images 

function [lesion_img, lesion_mask] = place_lesion(x,y, bg_img, tumor_shape_data, tumor_diameter, tumor_volume, heterogenous, max_iter, lesion, mask)

if nargin < 9
    [lesion, mask] = generate_tumor(tumor_shape_data, tumor_diameter, tumor_volume, heterogenous, max_iter);
end

f = figure(1);
% Place tumor
height = size(lesion, 1);
width = size(lesion,2);
x_ = round(x - width/2);
y_ = round(y - height/2);
lesion_img = zeros(128);
lesion_img(y_:height+y_-1, x_:width+x_-1) = lesion;
lesion_mask = zeros(128);
lesion_mask(y_:height+y_-1, x_:width+x_-1) = mask;

subplot(1,2,1); imagesc(bg_img+lesion_img); colorbar; colormap gray; 
if max(bg_img(:)) > 10
    caxis([0 10]);
end
axesHandles = findobj(get(f,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

subplot(1,2,2); imagesc(lesion_mask); colorbar; colormap gray; 
axesHandles = findobj(get(f,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

str = input('Keep tumor (y/n)?\n','s'); 
if ~strcmp(str,'y')
    [lesion_img, lesion_mask] = place_lesion(x,y, bg_img, tumor_shape_data, tumor_diameter, tumor_volume, heterogenous, max_iter); 
else
    str = input('Would you like to rotate the tumor (y/n)?\n','s');
    if strcmp(str,'y')
        for angle = 0:15:360
            [lesion_img,lesion_mask] = get_tumor_orientation(x, y, bg_img, lesion, mask, angle);
            str = input('Is this a good angle (y/n)?\n','s');
            if strcmp(str,'y')
                break; 
            end
        end
    end
    str = input('Adjust tumor placement (y/n)?\n','s');
    if strcmp(str,'y')
        [x,y] = get_tumor_position(bg_img);
        [row,col] = ind2sub(size(lesion_mask),find(lesion_mask>0));
        lesion = lesion_img(min(row):max(row), min(col):max(col));
        mask = lesion_mask(min(row):max(row), min(col):max(col));
        [lesion_img, lesion_mask] = place_lesion(x,y, bg_img, tumor_shape_data, tumor_diameter, tumor_volume, heterogenous, max_iter, lesion, mask);
    end
end

end

function [x,y] = get_tumor_position(bg_img)
f = figure(1);
imagesc(bg_img); colorbar; colormap gray;
if max(bg_img(:)) > 10
    caxis([0 10]);
end
axesHandles = findobj(get(f,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');
fprintf('Provide tumor location\n'); 
[x,y] = ginput(1);
x = round(x); y = round(y); 
end

function [lesion_img,lesion_mask] = get_tumor_orientation(x, y, bg_img, lesion, mask, angle) 

lesion_rot = imrotate(lesion, angle);
mask_rot = imrotate(mask, angle);

% Place tumor
height = size(lesion_rot, 1);
width = size(lesion_rot,2);
lesion_img = zeros(128);
lesion_img(y:height+y-1, x:width+x-1) = lesion_rot;
lesion_mask = zeros(128);
lesion_mask(y:height+y-1, x:width+x-1) = mask_rot;

f = figure(1);
subplot(1,2,1); imagesc(bg_img+lesion_img); colorbar; colormap gray; 
if max(bg_img(:)) > 10
    caxis([0 10]);
end
axesHandles = findobj(get(f,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

subplot(1,2,2); imagesc(lesion_mask); colorbar; colormap gray; 
axesHandles = findobj(get(f,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

end

%% Helper functions to generate random tumors from patient tumor statistics 

function [output, lesion_mask] = generate_tumor(tumor_shape_data, tumor_diameter, tumor_volume, heterogenous, max_iter)
% Generates random tumors of different shape, size, intensity and
% heterogenity based on patient data 

% Inputs: 
% tumor_shape_data - 4x5xPatients (using 5 harmonic coefficients) 
% tumor_size - patient tumor diameter data 
% tumor_intensity, tumor_variance - patient tumor intensity/variance data 
% heterogenous - boolean, if true then generate heterogenity
% lower intensity in the center of the lesion 
% max_iter - maximum allowed iterations to attain a resonable tumor shape 

% Outputs: 
% output - 
% lesion_mask -  binary mask of tumor 

% Generate random tumor shape 
counter = 0; 
while counter < max_iter % stop generating tumors after max_iter or after you get a valid tumor shape
    counter = counter + 1; 
    rFSDs = zeros(4,5); % generated Efourier coefficients
    for x = 1:4
        for y = 1:5
            tmp_seq = squeeze(tumor_shape_data(x,y,:));
            pd = fitdist(tmp_seq,'kernel'); 
            if all(tmp_seq==tmp_seq(1))
                rFSDs(x,y) = tmp_seq(1);
            else
                % Kernel dist
                rFSDs(x,y) = random(pd); 
            end
        end
    end
    outln = rEfourier(rFSDs, size(rFSDs,2), 100);
    outln = [outln; outln(1,:)];
    outln = abs(flipud(outln.*50)); % increase resolution
    bw = poly2mask(outln(:,2), outln(:,1), 41*50, 31*50);
    
    if sum(bw(:))==0
        disp('bw is zero');
    end
    
    % Check for self intersection
    [x0,~,~] = selfintersect(outln(:,2),outln(:,1));
    if isempty(x0)
        [I,J] = ind2sub(size(bw),find(bw>0));
        output = bw(min(I):max(I), min(J):max(J));
        break; 
    end
end

if counter == max_iter
    disp('ERROR: Too many iterations.');
    if ~exist('output')
        output = zeros(100); 
    end
end

% Create lesion mask 
lesion_mask = output > 0; 

% Model heterogenity in tumor shape 
if heterogenous == 1 
    output = generate_tumor_heterogeneity(output, 0.2, 0.8); 
end

% Apply random rotation 
angle = unifrnd(0,360);
output = imrotate(output, angle);
lesion_mask = imrotate(lesion_mask, angle);
[row,col] = ind2sub(size(lesion_mask),find(lesion_mask>0));
output = output(min(row):max(row), min(col):max(col));
lesion_mask = lesion_mask(min(row):max(row), min(col):max(col));
% fprintf('Angle of rotation: %.2f\n',angle);

% Model tumor size keeping the tumor volume within the bounds of the real
% patient data 
diameter = generate_tumor_diameter(tumor_diameter);
[x,y] = size(output);
scaling_factor = sqrt(diameter^2/(x^2+y^2));
x = ceil(scaling_factor*x);
y = ceil(scaling_factor*y);
lesion_mask_resize = imresize(lesion_mask, [x,y]);
while (sum(lesion_mask_resize(:)) < min(tumor_volume)) || (sum(lesion_mask_resize(:)) > max(tumor_volume))
    if sum(lesion_mask_resize(:)) < min(tumor_volume)
        disp('Tumor size too small');
    elseif sum(lesion_mask_resize(:)) > max(tumor_volume)
        disp('Tumor size too large');
    end
    
    diameter = generate_tumor_diameter(tumor_diameter);
    [x,y] = size(output);
    scaling_factor = sqrt(diameter^2/(x^2+y^2));
    x = ceil(scaling_factor*x);
    y = ceil(scaling_factor*y);
    lesion_mask_resize = imresize(lesion_mask, [x,y]);
end
lesion_mask = lesion_mask_resize; 
output = imresize(output, [x,y]); 
output = round(output); 
output = output .* lesion_mask; 

% Remove heterogenity if lesion mask is too small 
if sum(lesion_mask(:)) <= 100
    output(output==2) = 1; 
end

% Check if there are multiple objects 
cc = bwconncomp(lesion_mask); 
if cc.NumObjects ~= 1
    [output, lesion_mask] = generate_tumor(tumor_shape_data, tumor_diameter, tumor_volume, heterogenous, max_iter); 
    disp('Multiple ojects in tumor'); 
end

end

function output_mask = generate_tumor_heterogeneity(tumor_mask, size_low, size_high)
% Generate heterogeneity in tumor 
    tmp_large = tumor_mask;
    rnd_num = unifrnd(size_low,size_high); % shrink by randomly chosen percentage
    dim = round(size(tmp_large) .* rnd_num);
    tmp_small = imresize(tmp_large, dim);
    tmp_small = tmp_small>0;
    dim2 = size(tmp_small);
    dim = size(tmp_large);
    diff = round(abs(dim-dim2)./2);
    tmp_zero = zeros(size(tmp_large));
    tmp_zero(diff(1):diff(1)+dim2(1)-1, diff(2):diff(2)+dim2(2)-1) = tmp_small;
    output_mask = tmp_large + tmp_zero .* tmp_large; % get rid of overlap and sum to add heterogenous region
end

function output_mask = generate_tumor_intensity(tumor_mask, mean_tumor_intensity, tumor_variance, mean_background_intensity, background_variance)
% Generate variance in tumor intensity 

dim = find(tumor_mask>0);
[y, x] = ind2sub(size(tumor_mask), dim);
for i = 1:length(y)
    for j = 1:length(x)
        if tumor_mask(y(i),x(j)) == 1
            low = mean_tumor_intensity * (1-sqrt(mean(tumor_variance))/10); high = mean_tumor_intensity * (1+sqrt(mean(tumor_variance))/10);
            lesion_new_activity = max(min(mean([low, high]) + (high-low)/5 * randn(1), high), low);  % randomly change each pixel intensity in lesion with gaussian
            tumor_mask(y(i),x(j)) = lesion_new_activity;
        elseif tumor_mask(y(i),x(j)) == 2
            low = mean_background_intensity * (1-sqrt(mean(background_variance))/10); high = mean_background_intensity * (1+sqrt(mean(background_variance))/10);
            bg_new_activity = max(min(mean([low, high]) + (high-low)/5 * randn(1), high), low);  % randomly change each pixel intensity in lesion with gaussian
            tumor_mask(y(i),x(j)) = bg_new_activity; 
        end
    end
end
tumor_mask(tumor_mask<0) = 0; % constrain to positive values 
output_mask = tumor_mask;

end

function diameter = generate_tumor_diameter(tumor_diameter)
% Generate tumor diameter with kernel distribution 

pd = fitdist(tumor_diameter,'kernel');
diameter = random(pd); 

if diameter < min(tumor_diameter)
    diameter = generate_tumor_diameter(tumor_diameter); 
end

end

function TBR = generate_TBR(tumor_bg_intensity)
% Generate tumor to background ratio with kernel distribution 

pd = fitdist(tumor_bg_intensity,'kernel');
TBR = random(pd); 

if TBR < min(tumor_bg_intensity) % keep the TBR from being too low 
    TBR = generate_tumor_diameter(tumor_bg_intensity); 
end

end

