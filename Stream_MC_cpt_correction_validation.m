%% Monte Carlo Simulation and Validation with Stream Control Points

%% Import data
load('IR_stream_cpts_2016.mat')% IR temps at control point locations; blurred images removed (rev)
load('C:\Users\Emily\Google Drive\Peru_2016\Peru 2016 Matlab Data\Stream_cpts_2016.mat') %in Peru 2016 Matlab Data; ibutton temps

%% Create data vectors of equal length and timing
IR_times_2016 = IR_times_2016_rev(21:end,:);%get rid of IR data before ibutton data
IR_stream_cpts_2016 = IRstream_cpts_2016(:,21:end);

t1=datetime(2016,08,05,15,40,00); %starting time
t2=datetime(2016,08,09,08,30,00); %ending time
ibutton_time=(t1:hours(.25/3):t2)'; %datetime variable
[ib_time,idx] =intersect(ibutton_time,IR_times_2016); % find times where there is ibutton and IR data   

ibutton_stream_temps_2016 = Stream_cpts_2016(idx,:)'; %uses index to find ibutton temps at IR times

% load('IR_iB_2016_synchronous_data.mat') %saved IR and ibutton data from 2016 that overlaps

% plot IR and ibutton data to ensure they are the same length and times
plot(ib_time,mean(ibutton_stream_temps_2016,1),'b.'); hold on;
plot(IR_times_2016, mean(IR_stream_cpts_2016,1),'r.')
xlabel('Time');ylabel('Temperature (^oC)');title('Overlapping IR and ibutton Data')
legend('iButtons','IR Camera')
xlim([min(datenum(IR_times_2016)) max(datenum(IR_times_2016))])
%% Calculate Residuals
stream_res = ibutton_stream_temps_2016 - IR_stream_cpts_2016; % temperature residuals for stream

figure
plot(IR_times_2016,stream_res)
xlabel('Time');ylabel('Residual (^oC)');title('Stream Temperature Residuals 2016')
xlim([min(datenum(IR_times_2016)) max(datenum(IR_times_2016))])
%% Assignment of model parameters
n_stream = 1; % number of stream control poins   -value from 1 to 16

Boundaries = [1 13]; % Assign feasible ranges: n_stream
% need to somehow exlude combinations that choose no control points

%% UNIFORM RANDOM SAMPLING ROUTINE

NS = 1000; % Sample size  
[NP,NB]=size(Boundaries); % Check number of parameters NS

Cpt_quant = zeros(NS,NP); %how many of each type of control point will be used
for i=1:NP % Create random sample for each parameter (Monte Carlo); each column is one parameter
    Cpt_quant(:,i) = randi([Boundaries(i,1) Boundaries(i,2)],1,NS);
end

%% Run model for all random parameter sets
loc_index = cell(NS,1); % creates cell array to save row locations of chosen control points
stream_cp_temps = cell(NS,1); % creates cell array to save the stream temps of these chosen control points
valid_IR_Temps = cell(NS,1); % creates cell array to save the validation cpt temps
valid_locs = cell(NS,1); % creates cell array to save the validation cpt locations

count = 1; % Model Run #
for j=1:NS
    %% Model Discharge Inputs pulled from Monte Carlo random sample
    n_stream(j,1) = Cpt_quant(j); % number of stream control points
    
    % Randomly choose which rows (control point location residuals)
    [stream_locs, idx]= datasample(stream_res,n_stream(j),1,'Replace',false); % residuals from stream cpts
    loc_index{j} = idx; % saves locations of chosen control points
    stream_cp_temps{j} = stream_locs; % saves temps through time at these locations
    
    % Correct IR Temp using mean of residuals
    mean_corr_fact(:,:,j) = mean(stream_locs,1); % average the residuals from the chosen cpts to calculate correction factor
    for i = 1:14
    corr_Temp(i,:,j) = IR_stream_cpts_2016(i,:) + mean_corr_fact(:,:,j); % add the average residual to the IR temps (since IR temps are usually too cold)
    end
    
    % Calculate the root mean square error of locations not used for correction factor
    idx_range = (1:1:14); % possible range of index values
    valid_loc = idx_range(~ismember(idx_range,idx_range(idx))); % indexes not used for correction factor
    valid_locs{j} = valid_loc; % saved validation cpt locations
    
    valid_IR_temps = corr_Temp(valid_loc,:,j); % corrected IR temps from points not used for correction factor
    valid_IR_Temps{j} = valid_IR_temps; % saved temps at validation cpt location through time
    
    rmse(j,:)=sqrt(sum(sum((ibutton_stream_temps_2016(valid_loc,:)-valid_IR_temps).^2))/numel(ibutton_stream_temps_2016)); % rmse at of model at each time
    
    % Display model Run Completed
    X = [' Parameter Set Completed: ',num2str(count)]; disp(X)
    count = count+1;
end

%% Plot of corrected temperatures
% 
% figure % plots all average corrected IR temps from all simulations
% plot(IR_times_2016,mean(ibutton_stream_temps_2016,1),'k')
% hold on; plot(IR_times_2016,mean(IR_stream_cpts_2016,1),'r.')
% for run = 1:NS
% plot(IR_times_2016,mean(corr_Temp(:,:,run),1));hold on; 
% end
% hold on; plot(IR_times_2016,mean(ibutton_stream_temps_2016,1),'k','LineWidth',3)
% xlabel('Time'); ylabel('Temperature (^oC)'); title('Corrected Stream Temperatures from Using Different Stream Control Point Combinations')
% xlim([min(datenum(IR_times_2016)) max(datenum(IR_times_2016))])
% legend('In Stream Temp','IR Stream Temp')

%% Plots of residuals

figure % plot error vs number of control points
scatter(n_stream,rmse);  
xlabel('Number of Stream Control Points');ylabel('Mean RMSE')
title('RMSE Using Different Stream Control Point Combinations')

figure % plot error vs number of control points
boxplot(rmse,n_stream)
xlabel('Number of Stream Control Points');ylabel('Mean RMSE')
title('RMSE Using Different Stream Control Point Combinations')

%% Histogram of frequencies
figure % frequency of the number of control points used
histogram(n_stream)
xlabel('Number of Stream Control Points for Correction');ylabel('Frequency')