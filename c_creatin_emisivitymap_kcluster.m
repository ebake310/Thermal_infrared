%% Emissivity Map 
% By Caroline Aubry-Wake
% last edited January 17, 2015
% Image segmentation 
% This code runs aimage segmentation algorithm using k-means clustering to select the different type of
% surface. each region is then giben a certain emissivity, which is used in the planck function correction. 


%% Load Data

% Opening digital reference image 
load ('C:\Users\labusr\Documents\MATLAB\organized\DIG_registered.mat');

% Apply a sharpening filter to the image
I = imsharpen(DIG_registered);
figure, imshow(I)

%% Kmeans Clustering
% classify images regions using k-means clustering

% Apply a sharpening filter to the image
I = imsharpen(DIG_registered);
figure, imshow(I)

% classify images regions using k-means clustering
% this script is for intensity (grayscale) images. 
% It would have to be modified if color images were used instaed.  

ab = double(I(:,:,1)); % keep only a and b of lab colors
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,1); % gets all pixels of a anb b in a column 

nColors =5; % The number of colors that are wanted. 

% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors, 'distance', 'sqEuclidean');

% label pixels in the images from the results of KMeans  
% reput row and columns of ab in image setting                                
pixel_labels = reshape(cluster_idx,nrows,ncols);

% create images that segment imgae by colors                                                          
figure,imshow(pixel_labels,[]), title('image labeled by cluster index');                                  
  colormap (jet)
  colorbar

% Segment each color into a new images 
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 1]);

for p = 1:nColors
    color = I;
    color(rgb_label ~= p) = 0;
    segmented_images{p} = color;
end


 figure,imshow(segmented_images{1}), title('objects in cluster 1');
 figure,imshow(segmented_images{2}), title('objects in cluster 2');   
 figure,imshow(segmented_images{3}), title('objects in cluster 3');
 figure,imshow(segmented_images{4}), title('objects in cluster 4');
 figure,imshow(segmented_images{5}), title('objects in cluster 5');

% Visually analyze the results. Types can be combined to obtain the right
% cover. 
  
  % Here, 4 type are wanted:
  
  %% Defining areas
  % and give them emissivity values
em_map = zeros(768, 1024);
  
em_map(find(segmented_images{5})) = 0.98; % ice
em_map(find(segmented_images{1})) = 0.97; % snow
    em_map(find(segmented_images{4})) = 0.97; % snow
     em_map(find(segmented_images{3})) = 0.98; % broken ice
      em_map(find(segmented_images{2})) = 0.9;% rock
    imagesc(em_map); colormap(map); colorbar;


save('em_map', 'em_map');


%% Make figure for emissivity map %%
% saving emissivity map figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
fname = 'C:\Users\labusr\Documents\MATLAB\organized\C_Emissivity';
save(fullfile(fname, 'em_map'));
saveas(gca, fullfile(fname, 'em_map'), 'tif');


%% add glacier outline to figure
 load('C:\Users\labusr\Documents\MATLAB\organized\H_defining_regions\glaxcier_outline.mat')
% my glacier outline was already defined manually on the digital image
width = 7;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 8;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
lw_zero = 0.5;

axis1 = 262;
axis2 = 295;

axis3 = 262;
axis4 = 280;

 [~, threshold] = edge(outline1, 'sobel');
 fudgeFactor = 2;
 M = edge(outline,'sobel', threshold * fudgeFactor);

 imshow(M)
 M = double(M);
 se = strel('ball',3,3);
M2 = imdilate(M,se);
 imagesc(M2)
 
 A = roipoly;
 M2(A) = 0;
 A = find(M2>3.3);
 em_map2=em_map;
 em_map2(A) = 0;
imagesc(em_map2); colormap(map);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties
fname = 'C:\Users\labusr\Documents\MATLAB\organized\C_Emissivity';
saveas(gca, fullfile(fname, 'em_map2'), 'tif');


% -------------------------------------------------------------------------
% %% find sub element of a cluster 
% % Not needed in this case. Would be needed to differentiate light debris
% % from thick debris?
% 
% % select only light object in cluster using the L layer. 
% 
% % for each value, find the mean of the center value of row a and b 
% %and sort them (1 is rock, 2 is blue, 3 is white)
% mean_cluster_value = mean(cluster_center,2);
% [tmp, idx] = sort(mean_cluster_value);
% % for each clutser, extract light colorrs (white)
% a_cluster_num = idx(1);
% b_cluster_num = idx(2);
% c_cluster_num = idx(3);
% d_cluster_num = idx(4);
% e_cluster_num = idx(5);
% 
% % select only the high luminosity for each cluster
% 
% % show cluster A
% L = lab_I(:,:,1); % for the luminosity layer
% a_idx = find(pixel_labels == a_cluster_num); % all the blue pixels
% L_a = L(a_idx); % blue pixels in the luminosity layer
% is_light_a = im2bw(L_a,0.1); %% change threshodl value
% nuclei_labels = repmat(uint8(0),[nrows ncols]);
% nuclei_labels(a_idx(is_light_a==true)) = 1;
% nuclei_labels = repmat(nuclei_labels,[1 1 3]);
% a_nuclei = I;
% a_nuclei(nuclei_labels ~= 1) = 0;
% figure, imshow(a_nuclei), title('light pixel from cluster a');
% 
% % cluster b
% L = lab_I(:,:,1);
% b_idx = find(pixel_labels == b_cluster_num);
% L_b = L(b_idx);
% is_light_b = im2bw(L_b,0.1);
% b_labels = repmat(uint8(0),[nrows ncols]);
% b_labels(b_idx(is_light_b==true)) = 1;
% b_labels = repmat(b_labels,[1 1 3]);
% b_nuclei = I;
% b_nuclei(b_labels ~= 1) = 0;
% figure, imshow(b_nuclei ), title('light pixel from cluster b');
% 
% % cluster c
% L = lab_I(:,:,1);
% c_idx = find(pixel_labels == c_cluster_num);
% L_c = L(c_idx);
% is_light_c = im2bw(L_c,0.3);
% c_labels = repmat(uint8(0),[nrows ncols]);
% c_labels(c_idx(is_light_c==true)) = 1;
% c_labels = repmat(c_labels,[1 1 3]);
% c_nuclei = I;
% c_nuclei(c_labels ~= 1) = 0;
% figure, imshow(c_nuclei ), title('light pixel from cluster c');
% 
% % cluster d
% L = lab_I(:,:,1);
% d_idx = find(pixel_labels == d_cluster_num);
% L_d = L(d_idx);
% is_light_d = im2bw(L_d,0.25);
% d_labels = repmat(uint8(0),[nrows ncols]);
% d_labels(d_idx(is_light_d==true)) = 1;
% d_labels = repmat(d_labels,[1 1 3]);
% d_nuclei = I;
% d_nuclei(d_labels ~= 1) = 0;
% figure, imshow(d_nuclei), title('light pixel from dluster d');
% 
% %cluyster e
% L = lab_I(:,:,1);
% e_idx = find(pixel_labels == e_cluster_num);
% L_e = L(e_idx);
% is_light_e = im2bw(L_e,0.3);
% e_labels = repmat(uint8(0),[nrows ncols]);
% e_labels(e_idx(is_light_e==true)) = 1;
% e_labels = repmat(e_labels,[1 1 3]);
% e_nuclei = I;
% e_nuclei(e_labels ~= 1) = 0;
% figure, imshow(e_nuclei), title('light pixel from cluster e');
% 
% %% More transformation on regions; 
% % adding up each light clutser
% glacier = brown_nuclei + white_nuclei;
% S = im2bw(glacier,0.01);
% imshow(S)
% 
% P = 1000000;
% S= bwareaopen(S, P, 4);% witgh 4 connectiovity
% imshow(S);
% 
% % find and fill small holes
% filled = imfill(S, 'holes');
% holes = filled & ~S;
% bigholes = bwareaopen(holes, 100);
% imshow(bigholes);
% smallholes= holes & ~bigholes;
% imshow(smallholes)
% filled = S  | smallholes;
% imshow(filled)
% 
% 
% 
% [~, threshold] = edge(filled, 'sobel');
% fudgeFactor = 1.1;
% M = edge(filled,'sobel', threshold * fudgeFactor);
% figure, imshowpair(M,I), title('edge, glacier1');
% figure, imshow(I)
% contour_glacier1 = M;
% shape_glacier1=glacier ;
% 
% save('glacier1', 'shape_glacier1','contour_glacier1');
