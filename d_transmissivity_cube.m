
%% Create a transmissivity cube
% By Caroline Aubry-Wake
% Created Jan 22
% Revised by Emily Baker 7/14/2017

load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\D_Transmisivity\Transmissivity_data.mat')

%%
% This script create a data cube of transmissivity values to then be used in the image correction script. flor the thermacam software, I pull out the transmisivity values for -5C,
% 0 and 5C, for RH of 30 to 99, and for distance 800 to 1800. I organized
% those into arrays and create a transmissivity data cube

% I obtained those tr values from the FLIR software

temp_range = [-5 5 10 15 20 25]+273.15; % atmospheric temps in Kelvin
RH_range = [10 30 40 50 60 70 80 95];
distance_range = [1 10 50 100 200 300 400 500 600 700 800 900 1000]; % the data point I have data for

[X, Y, Z] = meshgrid (temp_range,RH_range, distance_range);   
% The actual data V is create by building a data cube with the temp in X,
% the RH in Y and the distance in Z. and the transmissivity data in the
% cube. I created it by copy pasting the data I have  obtained from flir.

% TR_values i obtained by creating a cuve ok transmissivity values
% extracted from the FLIR software. I copy paste each grid in the cube. 
D1000m=xlsread('FLIR_transmissivity_values.xlsx','1000m');
D1000=D1000m(2:7,3:10);
D900m=xlsread('FLIR_transmissivity_values.xlsx','900m');
D900=D900m(2:7,3:10);
D800m=xlsread('FLIR_transmissivity_values.xlsx','800m');
D800=D800m(2:7,3:10);
D700m=xlsread('FLIR_transmissivity_values.xlsx','700m');
D700=D700m(2:7,3:10);
D600m=xlsread('FLIR_transmissivity_values.xlsx','600m');
D600=D600m(2:7,3:10);
D500m=xlsread('FLIR_transmissivity_values.xlsx','500m');
D500=D500m(2:7,3:10);
D400m=xlsread('FLIR_transmissivity_values.xlsx','400m');
D400=D400m(2:7,3:10);
D300m=xlsread('FLIR_transmissivity_values.xlsx','300m');
D300=D300m(2:7,3:10);
D200m=xlsread('FLIR_transmissivity_values.xlsx','200m');
D200=D200m(2:7,3:10);
D100m=xlsread('FLIR_transmissivity_values.xlsx','100m');
D100=D100m(2:7,3:10);
D50m=xlsread('FLIR_transmissivity_values.xlsx','50m');
D50=D50m(2:7,3:10);
D10m=xlsread('FLIR_transmissivity_values.xlsx','10m');
D10=D10m(2:7,3:10);
D1m=xlsread('FLIR_transmissivity_values.xlsx','1m');
D1=D1m(2:7,3:10);

TR_values = cat(3,D1,D10,D50,D100,D200,D300,D400,D500,D600,D700,D800,D900,D1000);
TR_values2=permute(TR_values,[2 1 3]);
temp_query = [258.15:0.01:308.15];
RH_query = [10:1:100];
distance_query = [1:1:800];
[Xq,Yq,Zq] = meshgrid(temp_query,RH_query, distance_query); % the points I want to have data for
%TR_cube= griddata(X,Y,Z,TR_values2,Xq,Yq,Zq); %the final cube that can be queried
TR_cube= interp3(X,Y,Z,TR_values2,Xq,Yq,Zq); %the final cube that can be queried

save TR_cube temp_range RH_range distance_range TR_values TR_cube -v7.3

% see some results 

% figure
% slice(Xq,Yq,Zq,TR_cube,20,40,0);
% colormap hsv
% shading flat
% 
% figure
% pcolor( temp_query,RH_query, TR_cube(:,:, 164))
% shading flat
% pcolor(TR_cube(:,:, 164))

%% Distance grid  
rows_query = 1:768;
column_query = transpose(1:1024);
PixelDistance=xlsread('2016_pixel_distances.xlsx','Distances');
distance_googleearth=PixelDistance(:,3); %google earth extracted pixel distances
rr=round(PixelDistance(:,1)); %pixel row
cc=round(transpose(PixelDistance(:,2))); %pixel column
distance_grid = griddata(cc,rr,distance_googleearth,column_query,rows_query);

%saved as:
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\D_Transmisivity\2016_pixel_distance_grid.mat')

%% testing finding values %%%%%
% To make sure my cube is properly build
T_air_test = 273.15+15;
RH_test = 50;
distance_test = 500;

A = find (RH_query == RH_test)
B = find (abs(temp_query-T_air_test)<0.001)
C = find (distance_query == distance_test)

tr_test = TR_cube(A,B,C)%should equal 0.655826963 
