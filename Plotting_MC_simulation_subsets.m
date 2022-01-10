%% Plotting Results from MC Correction Simulations
load('/Users/emilybaker/Google Drive/IR Codes_Peru 2015/iButton_residual_corrections/Validation_up_to_16_cpts_10000_combos.mat')
% file containing resulting from 10000 different control point combinations

%% Plot n = 1

n1_corr = corr_Temp(:,:,find(Cpt_quant == 1)); % corrected temps using 1 cpt

figure % plots all average corrected IR temps from subset of simulations using n cpts
plot(redoTime,mean(subibStreamTemps,1),'k')
for i = 1:100
plot(redoTime,mean(n1_corr(:,:,i),1));hold on; 
end
hold on; plot(redoTime,mean(subibStreamTemps,1),'k','LineWidth',3)
hold on; plot(redoTime,mean(IRstreamtemps,1),'r.')
xlabel('Time'); ylabel('Temperature (^oC)'); title('Corrected Stream Temperatures Using One Stream Control Point')
xlim([min(datenum(redoTime)) max(datenum(redoTime))])
legend('In Stream Temp','IR Stream Temp')

%% Plot n = 4

n4_corr = corr_Temp(:,:,find(Cpt_quant == 4)); % corrected temps using 1 cpt

figure % plots all average corrected IR temps from subset of simulations using n cpts
plot(redoTime,mean(subibStreamTemps,1),'k')
for i = 1:15
plot(redoTime,mean(n4_corr(:,:,i),1));hold on; 
end
hold on; plot(redoTime,mean(subibStreamTemps,1),'k','LineWidth',3)
hold on; plot(redoTime,mean(IRstreamtemps,1),'r.')
xlabel('Time'); ylabel('Temperature (^oC)'); title('Corrected Stream Temperatures Using Four Stream Control Points')
xlim([min(datenum(redoTime)) max(datenum(redoTime))])
legend('In Stream Temp','IR Stream Temp')

%% Plot n = 8

n8_corr = corr_Temp(:,:,find(Cpt_quant == 8)); % corrected temps using 1 cpt

figure % plots all average corrected IR temps from subset of simulations using n cpts
plot(redoTime,mean(subibStreamTemps,1),'k')
for i = 1:15
plot(redoTime,mean(n8_corr(:,:,i),1));hold on; 
end
hold on; plot(redoTime,mean(subibStreamTemps,1),'k','LineWidth',3)
hold on; plot(redoTime,mean(IRstreamtemps,1),'r.')
xlabel('Time'); ylabel('Temperature (^oC)'); title('Corrected Stream Temperatures Using Eight Stream Control Points')
xlim([min(datenum(redoTime)) max(datenum(redoTime))])
legend('In Stream Temp','IR Stream Temp')


%% Plot n = 12

n12_corr = corr_Temp(:,:,find(Cpt_quant == 12)); % corrected temps using 1 cpt

figure % plots all average corrected IR temps from subset of simulations using n cpts
plot(redoTime,mean(subibStreamTemps,1),'k')
for i = 1:15
plot(redoTime,mean(n12_corr(:,:,i),1));hold on; 
end
hold on; plot(redoTime,mean(subibStreamTemps,1),'k','LineWidth',3)
hold on; plot(redoTime,mean(IRstreamtemps,1),'r.')
xlabel('Time'); ylabel('Temperature (^oC)'); title('Corrected Stream Temperatures Using Twelve Stream Control Points')
xlim([min(datenum(redoTime)) max(datenum(redoTime))])
legend('In Stream Temp','IR Stream Temp')