%% Extracting the whole 2016 stream  
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\B_registration\data_reg_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Aligned_and_Corrected_Images\IRdata_Aug2016_set1.mat', 'IR_times_set1')

%% Delineate Stream
% imagesc(data_reg_set1(:,:,1)); 
% hold on; plot(x2,y2,'k')
% hold on;plot(x1,y1,'k') %to check where the extracted line plots
% xt = CP_c+5; yt = CP_r;
% str = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
% text(xt,yt,str,'FontSize',8)

% [x1,y1,Raw_T1] = improfile %upper_stream_edge
% [x2,y2,Raw_T2] = improfile %lower_stream_edge
load('upper_stream_edge.mat')
load('lower_stream_edge.mat')

y1 = round(y1());
stream = data_reg_set1(:,:,1)-273.15;
for i = 1:1024
   [c index] = min(abs(i-x1));
   stream(1:y1(index),i) = NaN; 
end

y2 = round(y2);
for i = 1:1024
   [c index] = min(abs(i-x2));
   stream(y2(i):768,i) = NaN; 
end
figure; imagesc(stream)

for i = 452:457;
    stream(545:557,i) = NaN;
end
figure; imagesc(stream)

%% Choose stream pixels
% I = data_reg_set1(:,:,1)-273.15;
% figure; imagesc(I)
I_r = reshape(stream,[],768*1024);
eliminate = find(I_r>11.5);
for i = eliminate
I_r(1,i) = NaN;
end
I_cropped = reshape(I_r,768,1024);
figure; imagesc(I_cropped)
[nr,nc] = size(I_cropped);
pcolor([I_cropped nan(nr,1); nan(1,nc+1)]);
shading flat;  set(gca, 'ydir', 'reverse');

eliminate2 = isnan(I_r);
eliminate_all = find(eliminate2==1);
%% Extract and correct stream pixels for all images
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\iButton_stream_temps_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\IR_stream_cpts_2016_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\IR_Stream_Control_Points_16.mat')
load('C:\Users\Emily\Google Drive\Peru_2016\Peru 2016\iButtons_2016\iButton_calibration_correction_16.mat')

IR_rev = mean_IRpixeltemps(:,[23:452]);
iB_rev = subibStreamTemps+corr_16(1:13,:);
stream_res = IR_rev - iB_rev;

cp1 = stream_res(2,:);
cp2 = stream_res(7,:);
cp3 = stream_res(12,:);
cpt_mean_res = (cp1+cp2+cp3)/3;
for i=1:13
adj_IR_T(i,:) = IR_rev(i,:) - cpt_mean_res;
end
abs_error = abs(adj_IR_T - iB_rev);
mean_error = mean(reshape(abs_error,[1,13*430]));

I_stream = zeros(768,1024,430);
for k = 1:430
I2 = data_reg_set1(:,:,22+k)-273.15;
% figure; imagesc(I2)
I_r = reshape(I2,[],768*1024);
for i = eliminate_all
I_r(1,i) = NaN;
end
missing = find(I_r<-200);
for j = missing
I_r(1,j) = NaN;
end

I_r_corr = I_r - cpt_mean_res(1,k); % perform offset correction 

I_cropped = reshape(I_r_corr,768,1024);
% figure; imagesc(I_cropped)
I_stream(:,:,k) = I_cropped;
end

% corrected extracted stream saved as:
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\Extracted_corr_whole_stream.mat')

%% Make movie of whole stream 2016
% figure;imagesc(I_stream(:,:,100))
% [nr,nc] = size(I_stream(:,:,100));
% pcolor([I_stream(:,:,100) nan(nr,1); nan(1,nc+1)]);
% shading flat;  set(gca, 'ydir', 'reverse');
% c = colorbar; c.Label.String = 'Temperature (^oC)'; colormap(jet)
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\Extracted_corr_whole_stream.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\Res_adj_long_temp_line_set1.mat')

t=char(IR_Times);
v = VideoWriter('IR_stream_Video_2016'); %creates the file that you want to build the movie in
v.FrameRate = 4; %determines how fast each frame of the movie goes (x frames per second)
v.Quality = 100; %100 to 0 where 100 is highest quality
open(v)

for k = 1:430
    f=figure('position', [100, 100, 1024, 768]);
imagesc(I_stream(:,:,k))
[nr,nc] = size(I_stream(:,:,k));
pcolor([I_stream(:,:,k) nan(nr,1); nan(1,nc+1)]);
shading flat;  set(gca, 'ydir', 'reverse');
c = colorbar; c.Label.String = 'Temperature (^oC)'; colormap(jet)
caxis([median(reshape(I_stream(:,:,k),[],768*1024),'omitnan')-1.5 median(reshape(I_stream(:,:,k),[],768*1024),'omitnan')+1.5])

  str=t(k,:);
  dim=[.42 .66 .3 .3];
  annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','non','FontWeight','bold');
      drawnow
      pause(0.1)
      F=getframe(f); %captures the figure
  writeVideo(v,F);
  close all %closes the figure so that I don't run out of memory
end  
close(v);