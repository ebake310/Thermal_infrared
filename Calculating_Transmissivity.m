%% Transmissivity Equation from FLIR
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\D_Transmisivity\Transmissivity_data.mat')
%%
True_T = 36.8; % temp at 0.3 m reguarless of RH, Distance, and Air temp
True_W = 1*(5.67*10^-8)*(True_T+273.15)^4; % Radiant Energy in Wm-2K-4 for Stephan-Boltzmann Law

Temp_K = NewTempC+273.15; 
Measured_W = 1*(5.67*10^-8)*(Temp_K).^4; % Measured Radiant Energy

Trans = True_W./Measured_W; % Calculated transmissivity values

% saved values as Transmissivity Data

%% Plots
figure; plot(RH,Trans,'.'); xlabel('RH %');ylabel('Transmissivity');
figure; plot(Distance,Trans,'.'); xlabel('Distance (m)');ylabel('Transmissivity');
figure; plot(AtmTempC,Trans,'.'); xlabel('Atm Temp (^oC)');ylabel('Transmissivity');

figure
scatter3(RH,Distance,AtmTempC,40,Trans,'.')
xlabel('RH');ylabel('Distance');zlabel('Atm Temp')


