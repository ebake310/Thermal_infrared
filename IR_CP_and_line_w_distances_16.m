%% Extracting pixels from IR images along stream  
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\B_registration\data_reg_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Aligned_and_Corrected_Images\IRdata_Aug2016_set1.mat', 'IR_times_set1')

%% Manually Choose Points Along Stream to Extract IR temps from and save
% these points: c = column, r = row, Raw_T = pixel temp value
image(data_reg(:,:,12))
[tCP_c,tCP_r,tCP_Raw_T] = impixel(data_reg(:,:,12)); %pixel choose coordinates to extract from one image
CP_Raw_T = mean(CP_Raw_T,2); % gets rid of the three repeat columns
% Assign Control point names
CP_names = java_array('java.lang.String', 14);
CP_names(1) = java.lang.String('S1'); CP_names(2) = java.lang.String('S2'); CP_names(3) = java.lang.String('S3'); CP_names(4) = java.lang.String('S4');
CP_names(5) = java.lang.String('S5'); CP_names(6) = java.lang.String('S6'); CP_names(7) = java.lang.String('S7'); CP_names(8) = java.lang.String('S8');
CP_names(9) = java.lang.String('S9'); CP_names(10) = java.lang.String('S10'); CP_names(11) = java.lang.String('S11'); CP_names(12) = java.lang.String('S12');
CP_names(13) = java.lang.String('S13'); CP_names(14) = java.lang.String('S14');

%% Saved control points
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\IR_Stream_Control_Points_16.mat') % stream control point locations and names
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\Long_temp_line_set1.mat')

%% Extract Pixels along whole length of stream
figure; imagesc(data_reg_set1(:,:,1)); hold on; plot(CP_c,CP_r,'k.')
plot(x,y,'k') %to check where the extracted line plots
xt = CP_c+5; yt = CP_r;
str = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
text(xt,yt,str,'FontSize',8)

[x,y,Raw_T] = improfile

Raw_T = Raw_T-273.15;
figure;plot(x,y,'r-'); view([0 -90]); hold on; plot(CP_c,CP_r,'y.');  

%% Linear Interpolation of Distances
% first remove all points along longitudinal temp line outside of the
% control point coordinates

% % This doesn't work to interpolate:
% F1 = scatteredInterpolant(CP_c,CP_r,CP_Distances);
% Line_distances2 = F1(x,y);
% % Interpolate distances manually in excel for each segment
% 
% % Use to check difference between manual and function interpolation
% % Manual excel interpolation is better because distance is always increasing
% point = (1:1:970);
% plot(point,Line_distances,'b'); hold on; plot(point,Line_distances2,'r');

%% Extract the line of pixels from every image: use similar code as in VI reflection pixel extraction

% Extract Temp values from each registered IR image
for i = 1:452  
Temp_Line_Cube_set1(:,:,i) = impixel(data_reg_set1(:,:,i),x,y); %extract same pixels from other images; VI = name of visual image
end

Temp_Line_Cube_set1 = squeeze(Temp_Line_Cube_set1(:,1,:))-273.15; % data cube containing the temperatures along the whole reach through time
% rows are each location, column is temps, pages is image (time)
for i = 1:452
    for j = 1:851
if Temp_Line_Cube_set1(j,i)<-100
    Temp_Line_Cube_set1(j,i) = NaN
else Temp_Line_Cube_set1(j,i) = Temp_Line_Cube_set1(j,i);
end
    end
end 
% save as: load('Long_temp_line_set1.mat') % line of pixels along the whole stream
%% Plot raw longitudinal stream temperature through time and space
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\Long_temp_line_set1.mat')

% figure
% d=repmat(Line_distances_16,1,452);
% t=repmat(IR_times_set1,1,858)';
% s=surf(d,t,Temp_Line_Cube_set1)
% s.EdgeColor = 'none';
% xlabel('Distance Downstream (m)');ylabel('Time');zlabel('Raw Stream Temperature (^oC)')

figure; imagesc([1:1:452],Line_distances_16,Temp_Line_Cube_set1);ax = gca;
xlabel('Time');ylabel('Distance Downstream (m)')
c = colorbar; c.Label.String ='Raw Radiant Temperature (^oC)';colormap('jet')
ax.XTick = [13 49 85 120 156 191 215 251 287 322 357 393 427];
ax.XTickLabel = {'','0:00','','12:00','','0:00','','0:00','','12:00','','0:00',''};

%% Plot longitudinal stream temperature through time and space corrected using mean of 3 ibutton residuals
% load('/Users/emilybaker/Google Drive/IR Codes_Peru 2016/Longitudinal_analysis/Long_temp_line_set1.mat', 'Temp_Line_Cube_set1')
% load('/Users/emilybaker/Google Drive/IR Codes_Peru 2016/Extracting_control_points/iButton_stream_temps_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\Long_temp_line_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\iButton_stream_temps_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\IR_stream_cpts_2016_set1.mat')
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

for i=1:851
adj_long_IR_T_set1(i,:) = Temp_Line_Cube_set1(i,[23:452]) - cpt_mean_res;
end

load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\IR_Stream_Control_Points_16.mat')
figure; imagesc([1:1:430],Line_distances_16,adj_long_IR_T_set1);ax = gca;
hold on; plot(ones(13,1),CP_Distances(1:13,:),'k.')
xlabel('Time');ylabel('Distance Downstream (m)')
c = colorbar; c.Label.String ='Corrected Radiant Temperature (^oC)';colormap('jet')
ax.XTick = [13 49 85 120 156 191 215 251 287 322 357 393 427];%caxis([0 15])
ax.XTickLabel = {'','0:00','','12:00','','0:00','','0:00','','12:00','','0:00',''};
xt = ones(13,1)*5; yt = CP_Distances(1:13,:);
str = {'1','2','3','4','5','6','7','8','9','10','11','12','13'};
text(xt,yt,str,'FontSize',8)

%% Plot residuals of IR with distance downstream
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\iButton_stream_temps_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\IR_stream_cpts_2016_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\IR_Stream_Control_Points_16.mat')
load('C:\Users\Emily\Google Drive\Peru_2016\Peru 2016\iButtons_2016\iButton_calibration_correction_16.mat')

IR_rev = mean_IRpixeltemps(:,[23:452]);
iB_rev = subibStreamTemps+corr_16(1:13,:);
stream_res = IR_rev - iB_rev;
mean_res = mean(abs(stream_res),2);

figure; plot(CP_Distances(1:13,:),mean_res,'k.')
hold on; plot(CP_Distances(1:13,:),ones(13,1)*3.1,'kx');
      yt = (ones(13,1)*3.1)+0.02; xt = CP_Distances(1:13,:);
      str = {'1','2','3','4','5','6','7','8','9','10','11','12','13'};
      text(xt,yt,str,'FontSize',8)
ylabel('Residual (^oC)');xlabel('Distance Downstream (m)'); title('Uncorrected IR Residuals')

IR_rev = adj_IR_T;
iB_rev = subibStreamTemps+corr_16(1:13,:);
stream_res = IR_rev - iB_rev;
mean_res = mean(abs(stream_res),2);

figure; plot(CP_Distances(1:13,:),mean_res,'k.')
hold on; plot(CP_Distances(1:13,:),ones(13,1)*0,'kx');
      yt = (ones(13,1)*0)+0.02; xt = CP_Distances(1:13,:);
      str = {'1','2','3','4','5','6','7','8','9','10','11','12','13'};
      text(xt,yt,str,'FontSize',8)
ylabel('Residual (^oC)');xlabel('Distance Downstream (m)'); title('Corrected IR Residuals')

%% F9.Plot residuals of IR through time
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\iButton_stream_temps_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\IR_stream_cpts_2016_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Longitudinal_analysis\IR_Stream_Control_Points_16.mat')
load('C:\Users\Emily\Google Drive\Peru_2016\Peru 2016\iButtons_2016\iButton_calibration_correction_16.mat')

IR_rev = mean_IRpixeltemps(:,[23:452]);
iB_rev = subibStreamTemps+corr_16(1:13,:);
stream_res = IR_rev - iB_rev;
mean_res = mean(abs(stream_res),1);

figure; plot(IR_Times,mean_res,'r.')
ylabel('Residual (^oC)');xlabel('Time'); 

IR_rev = adj_IR_T;
iB_rev = subibStreamTemps+corr_16(1:13,:);
stream_res = IR_rev - iB_rev;
mean_res = mean(abs(stream_res),1);

hold on; plot(IR_Times,mean_res,'b.')
ylabel('Residual (^oC)');xlabel('Time'); 
legend('Uncorrected IR Residuals','Corrected IR Residuals')
xlim([min(IR_Times) max(IR_Times)]);

%% Standard deviation of stream temps through time
IR_std = std(adj_long_IR_T_set1);
iB_std = std(iB_rev);
figure;histogram(IR_std,'FaceColor','r','BinWidth',0.01); 
hold on; histogram(iB_std,'FaceColor','b','BinWidth',0.01)
legend('IR STD','iButton STD')
xlabel('Standard Deviation'); ylabel('Frequency')
