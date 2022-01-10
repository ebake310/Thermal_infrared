%% Monte Carlo Simulation: Using Different Control Points to Correct IR Temperatures

%% Import data
load('IR_iB_2016_synchronous_stream_data.mat') %stream IR and ibutton control points from 2016
load('iButton_ground_cpts_2016.mat')
load('IR_Ground_cpts_2016.mat')

iB_ground_cpts = [G1_2016,G3_2016,G4_2016];
IR_ground_cpts = [G1,G3,G4]';
%% Create data vectors of equal length and timing
IR_ground_times_2016 = IR_times_2016_rev(21:end,:);%get rid of IR data before ibutton data
IR_ground_cpts_2016 = IR_ground_cpts(:,21:end); %IR ground temps

t1=datetime(2016,08,05,15,40,00); %starting time
t2=datetime(2016,08,09,08,30,00); %ending time
ibutton_time=(t1:hours(.25/3):t2)'; %datetime variable
[ib_ground_time,idx] =intersect(ibutton_time,IR_ground_times_2016); % find times where there is ibutton and IR data   

ibutton_ground_temps_2016 = iB_ground_cpts(idx,:)'; %uses index to find ibutton temps at IR times

% load('IR_iB_2016_synchronous_stream_data.mat') %saved IR and ibutton data from 2016 stream temps that overlaps
% load('IR_iB_2016_synchronous_ground_data.mat') %saved IR and ibutton data from 2016 ground temps that overlaps

% plot IR and ibutton data to ensure they are the same length and times
plot(ib_ground_time,mean(ibutton_ground_temps_2016,1),'b.'); hold on;
plot(IR_ground_times_2016, mean(IR_ground_cpts_2016,1),'r.')
xlabel('Time');ylabel('Temperature (^oC)');title('Overlapping IR and ibutton Ground Temperature Data')
legend('iButtons','IR Camera')
xlim([min(datenum(IR_ground_times_2016)) max(datenum(IR_ground_times_2016))])

%% Calculate Residuals
subibStreamTemps = ibutton_stream_temps_2016;
IRstreamtemps = IR_stream_cpts_2016;
subpampacpts = ibutton_ground_temps_2016;
controlpoints = IR_ground_cpts_2016;

stream_res = subibStreamTemps - IRstreamtemps; % temperature residuals for stream
ground_res = subpampacpts - controlpoints; % temperature residuals for pampa locations

figure %plot residuals
plot(ib_time,mean(stream_res,1),'b'); hold on;
plot(ib_time,mean(ground_res,1),'g')
xlabel('Time');ylabel('Residual (^oC)');title('Mean Temperature Residuals 2016')
xlim([min(datenum(ib_time)) max(datenum(ib_time))])
legend('Stream Cpts','Ground Cpts')
%% Assignment of model parameters
n_ground = 1; % number of ground control points  -value from 0 to 10
n_stream = 1; % number of stream control poins   -value from 0 to 17

Boundaries = [0 3; 0 14]; % Assign feasible ranges: n_ground, n_stream
% need to somehow exlude combinations that choose no control points

%% UNIFORM RANDOM SAMPLING ROUTINE

NS = 1000; % Sample size  
[NP,NB]=size(Boundaries); % Check number of parameters NS

Cpt_quant = zeros(NS,NP); %how many of each type of control point will be used
for i=1:NP % Create random sample for each parameter (Monte Carlo); each column is one parameter
    Cpt_quant(:,i) = randi([Boundaries(i,1) Boundaries(i,2)],1,NS);
end

%% Run model for all random parameter sets
% Preallocate Arrays
loc_index_g = cell(NS,1); % creates cell array to save row locations of chosen control points
loc_index_s = cell(NS,1);
stream_cp_temps = cell(NS,1); % creates cell array to save the stream temps of these chosen control points
ground_cp_temps = cell(NS,1); % creates cell array to save the ground temps of these chosen control points

count = 1; % Model Run #
for j=1:NS
    %% Model Discharge Inputs pulled from Monte Carlo random sample
    n_ground(j,1) = Cpt_quant(j,1); % number of ground control points
    n_stream(j,1) = Cpt_quant(j,2); % number of stream control points
    
    % Randomly choose which rows (control point location residuals)
    [ground_locs, idx_g] = datasample(ground_res,n_ground(j,1),1,'Replace',false); % residuals from ground cpts
    loc_index_g{j} = idx_g; % saves locations of chosen control points
    
    [stream_locs, idx_s] = datasample(stream_res,n_stream(j,1),1,'Replace',false); % residuals from stream cpts
    loc_index_s{j} = idx_s; % saves locations of chosen control points
    
    combined_res = vertcat(ground_locs,stream_locs); % combine residuals from stream and ground cpts
    
    % Correct IR Temp using mean of residuals
    mean_corr_fact(j,:) = mean(combined_res,1); % average the residuals from the chosen cpts
    for i = 1:14
    corr_Temp(i,:,j) = IRstreamtemps(i,:) + mean_corr_fact(j,:); % add the average residual to the IR temps (since IR temps are usually too cold)
    end
    
    % Calculate the root mean square error of locations not used for correction factor
    s_idx_range = (1:1:14); % possible range of index values
    valid_loc = s_idx_range(~ismember(s_idx_range,s_idx_range(idx_s))); % indexes not used for correction factor
    valid_locs{j} = valid_loc; % saved validation cpt locations
    
    valid_IR_temps = corr_Temp(valid_loc,:,j); % corrected IR temps from points not used for correction factor
    valid_IR_Temps{j} = valid_IR_temps; % saved temps at validation cpt location through time
    
    rmse(j,:)=sqrt(sum(sum((subibStreamTemps(valid_loc,:)-valid_IR_temps).^2))/numel(subibStreamTemps)); % rmse at of model at each time
    
    % Display model Run Completed
    X = [' Parameter Set Completed: ',num2str(count)]; disp(X)
    count = count+1;
end
% plot(redoTime,ground_rows,'k');hold on; plot(redoTime,ground_res,'r') %to
% make sure its working
%% Plots
% 
% figure % plots all average corrected IR temps from all simulations
% plot(ib_time,mean(subibStreamTemps,1),'k')
% hold on; plot(ib_time,mean(IRstreamtemps,1),'r.')
% for run = 1:NS
% plot(ib_time,mean(corr_Temp(:,:,run),1));hold on; 
% end
% hold on; plot(ib_time,mean(subibStreamTemps,1),'k','LineWidth',3)
% xlabel('Time'); ylabel('Temperature (^oC)'); title('Corrected Stream Temperatures from Using Different Control Point Combinations')
% xlim([min(datenum(ib_time)) max(datenum(ib_time))])
% legend('In Stream Temp','IR Stream Temp')

%% Plots
figure % plot error vs number of control points
scatter(n_stream,n_ground,40,rmse,'filled');  
cb = colorbar; cb.Label.String = 'Mean RMSE'; colormap(jet)
xlabel('Number of Stream Control Points');ylabel('Number of Ground Control Points')

figure % plot error vs number of control points
scatter3(n_stream,n_ground,rmse,40,rmse,'filled'); 
xlabel('# Stream Cpts');ylabel('# of Ground Cpts')
zlabel('Mean RMSE')
cb = colorbar; cb.Label.String = 'RMSE'; colormap(jet)
