% Image registration, 2014 Infrared images
% By Caroline Aubry-Wake
% Created October 6, 2014
% Last edit: Jan  10, 2015
% ---------------
% to compensate for camera movement, infrared images are being
% registered against a reference image ( here the first infrared image of the series);
% Once images are all registered, their pixels should all match, allowing
% me to follow pixels thourghout the time lapse. 
% ------------------------

%Images should already be organized in a data cube
  
% load('IRdata_Aug2016_rev.mat'); % load the data cube

data = data_rev;
fixed = data(:,:,1);% the reference image is the first one of the series

% Attempt to modiy the IR image to be compared to visual image; not
% needed for now
% A = find(fixed < 200);
% fixed(A) = nan;
% ave = nanmean(nanmean(fixed));  
% ir_ave = fixed-ave;
% imagesc(ir_ave); colorbar
% ir_ave2 = (ir_ave*15)+150;imagesc(ir_ave2); colorbar
% A = isnan(ir_ave2);
% ir_ave2(A) = 0;
% ir_ave3 = ir_ave2/250;
% imshow(ir_ave3)
% fixed =ir_ave3;
% fixed2 = uint8(ir_ave2);imshow(fixed2);colorbar;


%%%%%%%%%%%%%%%% Digital image registration%%%%%%%%%%%%%%%%%%%%%%
%             if the reference (fixed)  image is not already in matlab format (if format is
%             double instead of uint8 or uint16), the image need to be reduced to a
%             intensity (black and white image) by picking only the last channel.
%   
%             fixed image can also be the digital image that will be used to outline
%             element of the glacier:
            
% to get the digital image registered
% moving_dig = imread ('C:\Users\labusr\Documents\MATLAB\organized\B_registration\AC062415.JPG'); % 
% moving2 = imresize (moving_dig, 0.5); % reducing the size (the resoltuion
% is too different between IR and digital image to be comparable, so i
% reduce the size first)
% imshow(moving2)
            % In that case, registration is only done once, to get that one
            % digital image to match the infrared ones. save manually at
            % the end, with name something like DIG_registered_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Registration
% Based on a example on matlab website
% set counter to keep track of how many images have been done
n = 1; 

for ii = 1:492; % the number of each image

    moving = data(:,:,ii);

 %% Registration
% 
% figure, imshowpair(moving, fixed, 'montage')
%  title('Unregistered')

% pick inital optimizationand metric
[optimizer,metric] = imregconfig('multimodal') ;
movingRegisteredDefault = imregister(moving, fixed, 'affine', optimizer, metric);
                % translation is only shifting the moving image, there is not scaling or
                % rotation. 
% 
% figure, imshowpair(movingRegisteredDefault, fixed)
%  title('A: Default registration');
% 
%Adjust initial radius
optimizer.InitialRadius = optimizer.InitialRadius/5;
movingRegisteredAdjustedInitialRadius= imregister(moving, fixed, 'affine', optimizer, metric);

%  figure, imshowpair(movingRegisteredAdjustedInitialRadius, fixed)
%   title('B:Adjusted InitialRadius')

% Increasing number of iterations
optimizer.MaximumIterations = 500;
movingRegisteredAdjustedInitialRadius500 = imregister(moving, fixed, 'affine', optimizer, metric);
% 
%   figure, imshowpair(movingRegisteredAdjustedInitialRadius500, fixed)
%   title('C: Adjusted InitialRadius, MaximumIterations = 300, Adjusted InitialRadius.')

% estimating and applying geometric transform
%tformSimilarity_JPG = imregtform(moving_JPG,fixed,'similarity',optimizer,metric);
tformSimilarity = imregtform(moving,fixed,'similarity',optimizer,metric);

Rfixed = imref2d(size(fixed));
movingRegisteredRigid = imwarp(moving,tformSimilarity,'OutputView',Rfixed);
   figure, imshowpair(movingRegisteredRigid, fixed);
   title('D: Registration based on similarity transformation model.');

   moving_registered = movingRegisteredRigid;
   
% figure, imshowpair(moving_registered, fixed, 'montage');
% title(['Registered Image  and IR' num2str(ii)]);

%% Saving the data in a new data cube
if n==1;
data_reg_2016 = zeros(768, 1024, 492);
else
end %concerned that it keeps remaking this zero matrix so tried to make it do it only at beginning yay fixed it

data_reg_2016(:,:, ii)= moving_registered;

 clear moving_registered_IR movingRegisteredAffineWithIC_IR...
delimiterIn headerlinesIn IR_image metric...
movingRegisteredDefault_IR optimizer.InitialRadius ...
movingRegisteredAdjustedInitialRadius_IR...
optimizer.MaximumIterations optimizer ...
movingRegisteredAdjustedInitialRadius500_IR tformSimilarity_IR ...
Rfixed movingRegisteredRigid_IR ...
moving_registered moving

n = n+1

close all
end

% saving the final data
save data_reg_2016 data_reg_2016 -v7.3