%% Temperature comparison
% By Caroline Aubry-Wake
% Created Oct 10; Edited Jan 17, 2015
%
% This code compared the results form the raw and corrected images. 
% It extract the surface temperature form the IR images and from the measured data at all control points. 
close all 
clear all

% Change directory to access files:
% cd (

%% Load Data
% Data cube, raw and modified
load('B_registration\data_compiled_bckp.mat')

% Registered digital image
load('B_registration\DIG_registered.mat')

% Measured data
load('Z_field_data\wx_glacier_14.mat', 'SurfaceT_gx14')
load('Z_field_data\Time_GxandIR.mat', 'Time_IR', 'Time_GX', 'Time_stake_edge', 'Time_Hobo')
load('Z_field_data\buttons_edge_14.mat', 'bottom6', 'bottom9')

load('Z_graph\colormap_black.mat')
load('F_Correlation_Wx\location.mat')

%% Select regions of interest
% select the points in the digital image where the control point are
% imagesc(data_mod5(:,:,200), [265 290]); colormap(map); colorbar; 
% wx_location =  roipoly(DIG_registered);(data_mod2(:,:,200));   % 540-530 , 3 large 
% stake6_location = roipoly(DIG_registered);    % 540-530 , 4 large A being the region of interest - the location of the weather station
% stake9_location = roipoly(DIG_registered);    % 540-530 , 5 large A being the region of interest - the location of the weather station
% 
% all_location = wx_location +stake6_location+stake9_location;
% imshowpair(DIG_registered, all_location ) % check to make sure its all good 
% 
% % save location
% save location wx_location stake6_location stake9_location

A = wx_location;
B = stake6_location;
C = stake9_location;

%% Extract data from IR images
% the for loop select the area A, B and C in each image, find the average and other
% statistics of the area, and then compiles it into a new matrix, with a
% column for each image

m=1
for ii =1:211% the number of my images in my files
%% using data cube
IR_image_mod =  data_compiled(:,:,ii);
IR_image_or =  data_compiled(:,:,ii);

% Changing nan to 0
D = isnan(IR_image_mod);
IR_image_mod(D) = 0;

E = isnan(IR_image_or);
IR_image_or(E) = 0;

% Temperature of the modified IR image
BA = A.* IR_image_mod;   % B being the image, with temperature only for the area selected (A)
 CA = find(BA);   % Giving the index (row, column) that have temperature values
BB = B.* IR_image_mod;   % B being the image, with temperature only for the area selected (A)
 CB = find(BB);   % Giving the index (row, column) that have temperature values
BC = C.* IR_image_mod;   % B being the image, with temperature only for the area selected (A)
 CC = find(BC);   % Giving the index (row, column) that have temperature values

 % Temperature of the raw IR images
FA = A.* IR_image_or;   % B being the image, with temperature only for the area selected (A)
 GA = find(FA);   % Giving the index (row, column) that have temperature values
FB = B.* IR_image_or;   % B being the image, with temperature only for the area selected (A)
 GB = find(FB);   % Giving the index (row, column) that have temperature values
FC = C.* IR_image_or;   % B being the image, with temperature only for the area selected (A)
 GC = find(FC);   % Giving the index (row, column) that have temperature values
 
 % Get average temperature oif each area
 S1_mod_A = nanmean(nanmean(BA(CA))); % finding the average temperature of that region
 S1_mod_B =  nanmean(nanmean(BB(CB))); % finding the average temperature of that region
 S1_mod_C =  nanmean(nanmean(BC(CC))); % finding the average temperature of that region
 S1_or_A = nanmean(nanmean(FA(GA))); % finding the average temperature of that region
 S1_or_B = nanmean(nanmean(FB(GB))); % finding the average temperature of that region
 S1_or_C = nanmean(nanmean(FC(GC))); % finding the average temperature of that region

 % Compilation
 T_mod_wx (:,ii) = S1_mod_A; 
 T_mod_stake6 (:,ii) = S1_mod_B; 
 T_mod_stake9 (:,ii) = S1_mod_C; 
T_or_wx(:,ii) = S1_or_A; 
T_or_stake6(:,ii) = S1_or_B; 
T_or_stake9(:,ii) = S1_or_C; 

m = m+1
end


% Put them all in one matrix
T_extract (:, 1) = T_mod_wx;
T_extract (:, 2) = T_mod_stake6;
T_extract (:, 3) = T_mod_stake9;
T_extract (:, 4) = T_or_wx;
T_extract (:, 5) = T_or_stake6;
T_extract (:, 6) = T_or_stake9;

%% Remove first night (cloudy)
% Some images are blurry and need to be removed (in interpolated in
% between?)
T_extract2(:,1:3) = T_extract(:, 1:3);
T_extract2(32:71, 1:3) = nan;

% T_extract(33:41, 1:3) = nan;
% T_extract(52:58, 1:3) = nan;
% T_extract(65:71, 1:3) = nan;
% T_extract(207:211, 1:3) = nan;
% 
%% Format field data
% to be able to do a scatter plot, it need to be a one on one format
T_measured(:, 1) = interp1(Time_GX, SurfaceT_gx14, Time_IR)+ 273.15;
T_measured(:, 2) = interp1(Time_stake_edge, bottom6, Time_IR)+ 273.15;
T_measured(:, 3) = interp1(Time_stake_edge, bottom9, Time_IR)+ 273.15;


% Decide offset correction to apply

offset_wx = nanmean(T_measured(:,1) - T_extract(:,1))
offset_st6 = nanmean(T_measured(:,2) - T_extract(:,2))
offset_st9 = nanmean(T_measured(:,3) - T_extract(:,3))
offset = (offset_wx+offset_st6+offset_st9) /3

offset2_wx = nanmean(T_measured(:,1) - T_extract2(:,1))
offset2_st6 = nanmean(T_measured(:,2) - T_extract2(:,2))
offset2_st9 = nanmean(T_measured(:,3) - T_extract2(:,3))
offset2 = (offset2_wx+offset2_st6+offset2_st9) /3


% apply offset
T_extract(:,1:3) = T_extract(:, 1:3) + offset;
T_extract2(:,1:3) = T_extract2(:, 1:3) + offset2;




% Cutting the field measurement to the right lenght, to cover same period
% as images
start_gx = find(Time_GX == Time_IR(1));
end_gx = find(Time_GX == Time_IR(end));
start_edge =  find(Time_stake_edge == Time_IR(1));
end_edge = find(Time_stake_edge == Time_IR(end));
start_hobo = find(Time_Hobo == Time_IR(1));
end_hobo = find(Time_Hobo == Time_IR(end));

T_measured_gx(1:length(start_gx:end_gx),1) = SurfaceT_gx14 (start_gx:end_gx, :)+273.15;
T_measured_edge(1:length(start_edge:end_edge),1) = bottom6 (start_edge:end_edge, :) +273.15;
T_measured_edge(1:length(start_edge:end_edge),2) =bottom9( start_edge:end_edge, :)+273.15;

RH_measured_GX (1:length(start_gx:end_gx),1) = RH_up_gx14 (start_gx:end_gx, :);
RH_measured_OBS (1:length(start_hobo:end_hobo),1) = OBS_RH (start_hobo:end_hobo, :);
Tair_measured_GX (1:length(start_gx:end_gx),1) = airT_up_gx14 (start_gx:end_gx, :)+273.15;
Tair_measured_OBS(1:length(start_hobo:end_hobo),1) = OBS_T (start_hobo:end_hobo, :);
%% Scatter plots

% Creatiing scatter plot to see thre match between measured and corrected
% infrared images
width = 4;     % Width in inches
height = 12;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 12;      % Fontsize
lw = 1.5;      % LineWidth
msz = 40;       % MarkerSize
lw_zero = 0.5;

xtra = 0.02

%% overall r value
Textract_all(1:211, 1) = T_extract2(:, 1);
Textract_all(212:211*2, 1) = T_extract2(:, 2);
Textract_all(211*2+1:211*3, 1) = T_extract2(:, 3);

Tmeasured_all(1:211, 1) = T_measured(:, 1);
Tmeasured_all(212:211*2, 1) = T_measured(:, 2);
Tmeasured_all(211*2+1:211*3, 1) = T_measured(:, 3);



r_all = corrcoef (Tmeasured_all(:, 1),Textract_all(:, 1), 'rows', 'pairwise')
scatter (Tmeasured_all, Textract_all)
% For wx
h2 = figure(1);

h2_1 = subplot(3,1,1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*200]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties

% create one-one line
xmin = 265; xmax = 277;
ymin = xmin; ymax = ymin;
clear one_one_line
one_one_line(:,1) = [xmin:xmax]; 
one_one_line (:,2) =  [xmin:xmax]; 

% plot scatter
plot(one_one_line(:,1),one_one_line(:,2), 'k');
hold on
scatter (T_measured(:, 1), T_extract(:, 1),msz, 'r', 'filled', 'diamond');
scatter (T_measured(:, 1), T_extract2(:, 1),msz, 'k', 'filled',  'diamond');


% Graph set-up
axis ([xmin xmax xmin xmax]);
r_wx = corrcoef (T_measured(:, 1),T_extract(:, 1), 'rows', 'pairwise')
r_wx_2 = corrcoef (T_measured(:, 1),T_extract2(:, 1), 'rows', 'pairwise')

ylabel ('Corrected', 'FontSize',fsz, 'FontWeight', 'bold');
xlabel ('Measured','FontSize', fsz, 'FontWeight', 'bold');
%title ('Wx', 'FontSize', fsz, 'FontWeight', 'bold');

axis([265 276 265 276])
% set(gca,'YTick',265:2:276); %<- Still need to manually specific tick marks
% set(gca,'XTick',266:2:276); %<- Still need to manually specific tick marks


set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties

% 
% fname = 'F:\organized\organized\F_Correlation_Wx';
% saveas(gca, fullfile(fname, 'scatter_wx'), 'tif');
% print(h2, fullfile(fname, 'scatter_wx'),'-depsc','-r300');

% savefig ('scatter_WX.fig');
% saveas (gcf,'scatter_WX.jpg'); 

%% stake6 

h2_2 = subplot(3,1,2);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*200]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties

% Creating one-one line
xmin = 265; 
xmax = 295;
ymin = xmin;
ymax =xmax;
clear one_one_line
one_one_line(:,1) = [xmin:xmax]; 
one_one_line (:,2) = [xmin:xmax]; 

% plot Scatter
plot(one_one_line(:,1),one_one_line(:,2), 'k');
hold on
scatter (T_measured(:, 2),T_extract(:, 2), msz,  'r', 'filled','diamond');
scatter (T_measured(:, 2),T_extract2(:, 2), msz,  'k', 'filled', 'diamond');

% Graph set-up
axis ([xmin xmax xmin xmax])
ylabel ('Corrected', 'FontSize', fsz, 'FontWeight', 'bold');
xlabel ('Measured','FontSize', fsz, 'FontWeight', 'bold');
r_stake6 = corrcoef (T_measured(:, 2),T_extract(:, 2),'rows', 'pairwise')
r_stake6_2 = corrcoef (T_measured(:, 2),T_extract2(:, 2),'rows', 'pairwise')

%title ('Stake6', 'FontSize', fsz, 'FontWeight', 'bold');

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties

% fname = 'E:\organized\organized\F_Correlation_Wx';
% saveas(gca, fullfile(fname, 'scatter_st6'), 'tif');
% print(h2, fullfile(fname, 'scatter_st6'),'-depsc','-r300');


%% stake 9
h2_3 = subplot(3,1,3);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*200]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties

% One-one line
xmin = 265; 
xmax = 300;
ymin = xmin;
ymax =xmax;
clear one_one_line
one_one_line(:,1) = [xmin:xmax]; 
one_one_line (:,2) = [xmin:xmax]; 

plot(one_one_line(:,1),one_one_line(:,2), 'k');
hold on

% Scatter plot
scatter (T_measured(:, 3),T_extract2(:, 3), msz, 'r', 'filled','diamond');
scatter (T_measured(:, 3),T_extract(:, 3),msz, 'k', 'filled', 'diamond'); 

% Graph set-up
axis ([xmin xmax xmin xmax])
ylabel ('Corrected', 'FontSize', fsz, 'FontWeight', 'bold');
xlabel ('Measured','FontSize', fsz, 'FontWeight', 'bold');
%title ('Stake9', 'FontSize', fsz, 'FontWeight', 'bold');

r_stake9 = corrcoef (T_measured(:, 3),T_extract(:, 3),'rows', 'pairwise')
r_stake9_2 = corrcoef (T_measured(:, 3),T_extract2(:, 3),'rows', 'pairwise')

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties


pos1 = get(h2_1, 'Position') % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 xtra]
    set(h2_1, 'Position',new_pos1 ) % set new position of current sub - plot
pos1 = get(h2_2, 'Position') % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 xtra]
    set(h2_2, 'Position',new_pos1 ) % set new position of current sub - plot
pos1 = get(h2_3, 'Position') % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 xtra]
    set(h2_3, 'Position',new_pos1 ) % set new position of current sub - plot
    
    fname = 'F:\organized\organized\F_Correlation_Wx';
saveas(gca, fullfile(fname, 'scatter_all'), 'tif');
print(h2, fullfile(fname, 'scatter_all'),'-depsc','-r300');

%% Time series plot
width = 5;     % Width in inches
height = 10;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 6;      % Fontsize
lw = 1;      % LineWidth
msz = 8;       % MarkerSize
lw_zero = 0.5;
xtra = 0.04;
% Wx
h2 = figure;

h2_1 = subplot(5,1,1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*200]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties

plot(Time_GX(start_gx:end_gx), smooth(Tair_measured_GX, 10), 'r');
hold on
plot(Time_Hobo (start_hobo:end_hobo), smooth(Tair_measured_OBS+273.15, 10), 'b');
plot(Time_IR, ones(1, 211)*273.15, 'k', 'LineWidth', lw_zero);
set(gca,'XTick',[Time_IR(1):0.1250:Time_IR(end)]);
datetick('x', '', 'keepticks')
axis([Time_IR(1) Time_IR(207) 270 286])
set(gca,'YTick',270:5:285); %<- Still need to manually specific tick marks


 h2_2 = subplot(5,1,2);

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*150]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties

plot(Time_GX(start_gx:end_gx), smooth(RH_measured_GX, 10), 'r');
hold on
plot(Time_Hobo (start_hobo:end_hobo),smooth( RH_measured_OBS, 10), 'b');
set(gca,'XTick',[Time_IR(1):0.1250:Time_IR(end)]);
datetick('x', '', 'keepticks')
axis([Time_IR(1) Time_IR(207) 38 82])
set(gca,'YTick',45:10:80); %<- Still need to manually specific tick marks

h2_3 = subplot(5,1,3);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*200]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties
grid on

plot(Time_GX(start_gx:end_gx), T_measured_gx(:,1), 'k');
hold on
plot(Time_IR, T_extract(:,4), 'b');
plot(Time_IR, T_extract(:,1),'--r');
plot(Time_IR, T_extract2(:,1),'r');
plot(Time_IR, ones(1, 211)*273.15,'k' , 'LineWidth', lw_zero);
set(gca,'XTick',[Time_IR(1):0.1250:Time_IR(end)]);
datetick('x', '', 'keepticks')
axis([Time_IR(1) Time_IR(207) 258 278])
set(gca,'YTick',260:5:275); %<- Still need to manually specific tick marks


h2_4 = subplot(5,1,4);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*200]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties
grid on
plot(Time_stake_edge(start_edge:end_edge), T_measured_edge(:,1),'k');
hold on
plot(Time_IR, T_extract(:,5),'b');
plot(Time_IR, T_extract(:,2),'--r');
plot(Time_IR, T_extract2(:,2),'r');
plot(Time_IR, ones(1, 211)*273.15, 'k', 'LineWidth', lw_zero);
set(gca,'XTick',[Time_IR(1):0.1250:Time_IR(end)]);
datetick('x', '', 'keepticks')
axis([Time_IR(1) Time_IR(207) 262 293])
set(gca,'YTick',265:5:290); %<- Still need to manually specific tick marks

h2_5 = subplot(5,1,5);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*200]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', lw); %<- Set properties

plot(Time_stake_edge(start_edge:end_edge), T_measured_edge(:,2),'k');
hold on
plot(Time_IR, T_extract(:,6), 'b');
plot(Time_IR, T_extract(:,3),'--r');
plot(Time_IR, T_extract2(:,3),'r');
plot(Time_IR, ones(1, 211)*273.15, 'k', 'LineWidth', lw_zero);
set(gca,'XTick',[Time_IR(1):0.1250:Time_IR(end)]);
datetick('x', 'HHAM', 'keepticks')
axis([Time_IR(1) Time_IR(207) 262 302])
set(gca,'YTick',265:10:300); %<- Still need to manually specific tick marks

pos1 = get(h2_1, 'Position') % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 xtra]
    set(h2_1, 'Position',new_pos1 ) % set new position of current sub - plot
pos1 = get(h2_2, 'Position') % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 xtra]
    set(h2_2, 'Position',new_pos1 ) % set new position of current sub - plot
pos1 = get(h2_3, 'Position') % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 xtra]
    set(h2_3, 'Position',new_pos1 ) % set new position of current sub - plot
pos1 = get(h2_4, 'Position') % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 xtra]
    set(h2_4, 'Position',new_pos1 ) % set new position of current sub - plot
pos1 = get(h2_5, 'Position') % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 xtra]
    set(h2_5, 'Position',new_pos1 ) % set new position of current sub - plot
    
    
    
fname = 'E:\organized\organized\F_Correlation_Wx';
saveas(gca, fullfile(fname, 'timeseries'), 'tif');
print(h2, fullfile(fname, 'timeseries'),'-depsc','-r450');



%% Creating Time axis


% % Create a timeline for the images 
% 
% Time_string = num2str(Time_IR_yyyyddhhmm);
% Time_str2 = strcat(Time_string(:, 5:6),'/06/', Time_string(:, 1:4), ', ',Time_string(:, 7:8), ':',Time_string(:, 9:10));
% Time_str3 = [Time_string(:, 5:6),'/06/', Time_string(:, 1:4), ',  ',Time_string(:, 7:8), ':',Time_string(:, 9:10)]
% Timeformat = 'dd/mm/yyyy,HH:MM';
% Time_datenum = datenum(Time_str2, Timeformat);
% Time_IR= Time_datenum;
% 
% 
% % Creating time array form Gx
% 
% Time_string = num2str(Time_Gx(2:24099));
% Time_str2 = strcat(Time_string(:, 1:2),'/',Time_string(:, 3:4),'/', Time_string(:, 5:8), ', ',Time_string(:, 9:10), ':',Time_string(:,11:12));
% Timeformat = 'dd/mm/yyyy,HH:MM';
% Time_datenum = datenum(Time_str2, Timeformat);
% 
% Time_GX=  Time_datenum;


