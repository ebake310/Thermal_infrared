 %% Image correction
% By Caroline Aubry-Wake and Emily Baker
% Last edit Jan 21, 2016 

% This code extract the surface temperature of each pixel of the image.
% It start with the registered images and applies the planck function to 
% each pixel to obtain the surface temperature.

% Since we don't have camera weather for 2016, we only use the estimated
% reflected temperatures to correct the data

% Data used:
% Emissivity 
% Transmsiivity
% Registered images in a data cube
% Measured temperature for environment and for atmopshere

%% loaddata 
% data cube
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\B_registration\data_reg_set1.mat');
%imagesc(data_reg_set1(:,:,2)) %registered set 1 data in kelvin
data_cube = data_reg_set1(:,:,17:end);

% emissivity grid
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2015\C_Emissivity\em_map.mat'); % 0.96 for water

% transmissivity cube
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\D_Transmisivity\TR_cube.mat'); % new grid using FLIR values

%distance grid
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\D_Transmisivity\2016_pixel_distance_grid.mat'); %created in D_transmissivity

% Field measurement for atm and environement
% Station_Weather=xlsread('Weather_Station_Data.xlsx','IR_Weather');
% Station_Temp=Station_Weather(:,2); %temp at weather station in degrees Celcius
% Station_RH=Station_Weather(:,3); %RH at weather station in percent
% %load('C:\Users\labusr\Documents\MATLAB\organized\Z_field_data\wx_glacier_14.mat', 'airT_up_gx14', 'RH_up_gx14') %weather station temp
% %load('C:\Users\labusr\Documents\MATLAB\organized\Z_field_data\hobo_14.mat', 'OBS_T') %camera location temp
% Camera_Weather=xlsread('Lascar_IRcam_Quilcay15_TempRH.xlsx','Abbr_Camera_Weather');
% Camera_Temp=Camera_Weather(:,1); %temp at camera in degrees Celcius
load('C:\Users\Emily\Google Drive\Peru_2016\Peru 2016 Matlab Data\Weather_data_2016.mat', 'RH', 'Temp','Weath_times')

% time steps
% Times=xlsread('iButton_stream_temps.xlsx','IR_Times');
% DateTime=Times(:,4);
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Aligned_and_Corrected_Images\IRdata_Aug2016_set1.mat', 'IR_times_set1'); %to get the datetime variable I imported from the matlab spreadsheet by clicking into the excel sheet in Matlab and changing the column to datetime then importing that column as a variable

% digitsl image
 %load('DIG_registered_1.mat')
% radiant power to temperature parameter fit % needed the first time you
% deal with a dataset, then you can just load  them
%load('C:\Users\labusr\Documents\MATLAB\Planck\polyfit_variable.mat')

% load reflected temperatures estimated using modified caroline's code
% load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\refl_T_array_16_K_88.mat')
% load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\reflected_temps_2016.mat', 'times')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Other\reflected_temps_C_e_of_0.7.mat', 'times')

%% Constants
% defining constant to use in the image correction 
h = 6.626*10^-34;               % planck constnat [Js]
c = 2.998*10^8;                 % speed of light  in vacuum [m/s]
k = 1.3806488*10^-23;           % Bolztman constant
Trefl_default = 20 + 273.15;    % temperature of reflected environment [K]
Tatm_default = 20 + 273.15;     % Temperature of the atmopshere (path) [K]
em_default = 1;                 % emissivity. In varicapture, default is 1.
tr_default = 1 ;                % transmissivity
wv = [7.5:0.1:14]*10^-6;        % begin to end of wavelenght capture by camera [um]

constant = (2* pi * h * c^2) ;  % to make calculation simpler

% Data cube to save files
data_mod = zeros(768, 1024, 436);

%% Field data
% the environement temperature used is from the camera location (OBS)
% the atmosphere temperature used is the air temperature measured at the
% glacier weather station

% The measured data is processed to fit in thetime step of the images.
% make arrays the same time steps as the IR data
A = ismember(Weath_times,IR_times_set1);
B = find(A==1);
times_2016 = Weath_times(B,1); %times data overlap
RH_2016 = RH(B,1); % ibutton stream temps
Air_T_2016 = Temp(B,1);

A = ismember(times_2016,times);
B = find(A==1);
coinc_times_2016 = times_2016(B,1); %times data overlap
RH_2016_2 = RH_2016(B,1); % ibutton stream temps
Air_T_2016_2 = Air_T_2016(B,1);
% Est_refl_T = refl_T_array_16;
Est_refl_T = Air_T_2016_2; % make reflected T air T

% convert from Celcius to Kelvin
Tatm_actual=round(Air_T_2016_2+273.15,2); % air temp
Trefl=Est_refl_T+273.15; % environment/reflected temp
RH_actual=round(RH_2016_2);

%% Reshaping emissivity and transmissivity
% transform the emissivity and transmissity grid into array

% find the nan in a distance grid 
distance_nan = isnan(distance_grid);
% em_image = em_map-0.08; % emissivity of 0.88 for camera angle of 70 degrees
em_image = em_map; % emissivity of 0.96 
em_image(distance_nan) = nan;

% rehsape distance and em

size_image = size (distance_grid);

em_array = reshape (em_image, 1, size_image(1)*size_image(2));
distance_array = reshape(distance_grid, 1, size_image(1)*size_image(2));
distance_array = round(distance_array);
%% Correlation between M and T
% will be needed later

T_all = [260:1:310]; % range of temperature expected by your object/images; if you make this range too big the function goes crazy

for iii = 1:length (T_all);
  T = T_all(iii);
  M= (constant./(wv.^5) * 1 ./ (exp( (h*c)./(wv*k*T)) - 1)) ;
  T_sum (iii)= trapz(wv,M);
end 

% i fit a curve to the data set
p = polyfit (T_sum,T_all, 2);
y = polyval(p, T_sum);
R = corrcoef(T_sum, y);
scatter (T_sum, T_all,'b'); %the actual relationship between T and M
hold on
scatter(T_sum, y, 'r'); %the fit relationship between T and M
legend('Actual','Fit','location','southeast')
xlabel('Radiance')
ylabel('Temperature (K)')
 
 %% Correct images     

tic % to know how long each image takes to process

for ii = 1:436; % the number of my images in my files (times of IR data)
    
    % open one image at a time, reshape as an array
    IR_image = data_cube(:,:,ii);
    IR_image(distance_nan) = nan;
    
    IR_array = reshape(IR_image, 1, size_image(1)*size_image(2)); % changing the image into an array
    
    Trefl_actual = Trefl;
    
    for z = 1:length(IR_array);
        
        temp_query = [258.15:0.01:308.15];
        RH_query = [10:1:100];
        distance_query = [1:1:800];
        
        
        A = find (RH_query==RH_actual(ii)); 
        B = find(abs(temp_query-Tatm_actual(ii))<0.001); %was not working
        C = find (distance_query == distance_array(z));
        
        tr_array = TR_cube(A, B, C);
        checkfile = isempty(tr_array);
        if checkfile == 1;
            tr_array = nan; %this is giving me Nan seemingly every time, not querying temp all the time
        end
        % Going through each pixel of the array and processing the pixel to obtain the object
        % temperature
        
        % the temperature of the object given by the IR image
        Tobj_raw = IR_array(z);
        
        % Obtain the raidant power as given by the IR image
        % (M is the total radiant power of a blackbody of given temperature T in
        % wavelenght interval (wv, wv+dy))
        M_obj= (constant ./ (wv.^5) .* 1 ./ (exp( (h*c)./(wv*k*Tobj_raw)) - 1)) ;
        
        
        % To find the area under the curve (integral), a numerical trapezoid method
        % is used. The integral is evaluated at each value of wv (every 0.1um)
        sum_Mobj(z) = trapz(wv, M_obj);
        
        % find radiance I mfor the measured
        I_meas (z)= sum_Mobj(z)./ pi;
        
        % (as default Tr and Em are 1, the measured and object radiance are the
        % same.)
        
        % From Imeas, find T_obj
        
        % Find radiant power M for each component using measured temperature
        M_reflected_mod = (constant ./ (wv.^5) .* 1 ./ (exp( (h*c)./(wv*k*Trefl_actual(ii,1))) - 1)) ;
        M_atm_mod= (constant ./ (wv.^5) * 1 ./ (exp( (h*c)./(wv*k*Tatm_actual(ii))) - 1));
        
        sum_Mrefl_mod(z) = trapz(wv,M_reflected_mod); %integrating to solve for M reflected
        sum_Matm_mod(z) =  trapz(wv,M_atm_mod); %integrating to solve for M of atmosphere
        
        I_refl_mod(z)= sum_Mrefl_mod(z)/ pi; %calculating reflected I
        I_atm_mod(z) = sum_Matm_mod(z)/ pi; %calculating atmospheric I
        
        I_obj_mod(z) = (I_meas(z) - (I_refl_mod (z)* (1-em_array(z)) * tr_array)  -  (I_atm_mod (z) * (1-tr_array)))/ (tr_array * em_array(z)); % area under the curve
        sum_Mobj_mod(z) = I_obj_mod(z)*pi;
        
    end

% I obtain a radiant power for each pixel of the image for the spectral
% range captured by the camera.

% I need to transfer this radiant power to temperature. I use the
% correlation between M and T. 
 
for x =  1:length(sum_Mobj_mod); 
    T_obj_mod(x) = polyval(p, sum_Mobj_mod(x));
end 

% reshapen image
T_obj_grid = zeros(768, 1024);
T_obj_grid = reshape(T_obj_mod, 768, 1024);
% change the name, save as a different filename


data_mod (:,:,ii) = T_obj_grid; 
ii
toc
end 

save data_mod_trad data_mod coinc_times_2016 -v7.3
% reflected temps = air temps; emissivity of 0.96

% save data_mod_revE data_mod coinc_times_2016 -v7.3
% reflected temps = air temps; emissivity of 0.88 based on 70 deg angle

% save data_mod_est88 data_mod coinc_times_2016 -v7.3
% data_mod_est88 was made using e = .88 (70 deg, 2016) and trans from FLIR, and estimated
% reflected temp by using carolines method to back calculated reflected T

% %histogram(data_mod)
% bad_values = find(data_mod< -1000);
%   
% a = min(data_mod);
% b = min(a); %looking for bad values
% %% Checking range of corrected data
% 
% data_mod_reshape=reshape(data_mod,[1 556007424]);
% data_max=max(data_mod_reshape);
% data_min=min(data_mod_reshape);  
