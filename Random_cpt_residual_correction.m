%% Correction using mean adj_IR_Residuals from a few control points
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\IR_stream_cpts_2016_set1.mat')
load('C:\Users\Emily\Google Drive\Peru_2016\Peru 2016 Matlab Data\Stream_cpts_2016.mat')

%% Make vectors same lengths of time
stream_cpt_times = stream_cp_times(14:2:1052);  % take ibutton data every 10 minutes
iB_Stream_2016 = Stream_cpts_2016(14:2:1052,:); 

IRstreamtemps = mean_IRpixeltemps(:,23:end); %shorten IR data to length of time we have iButton data for
IR_Times = IR_times_set1(23:end,:);

[C,ia] = intersect(stream_cpt_times,IR_Times,'rows'); %finds matching times

subibStreamTemps = iB_Stream_2016(ia,1:13)'; % stream ibutton temperatures at the times we have IR data for
figure; plot(stream_cpt_times,iB_Stream_2016,'k.'); hold on; plot(IR_Times,subibStreamTemps,'.'); %check that data were extracted correctly

figure; p1= plot(stream_cpt_times,iB_Stream_2016,'b.'); hold on
p2 = plot(IR_Times,IRstreamtemps,'r.')
xlabel('Time');ylabel('Temperature (^oC)');
xlim([min(IR_Times) max(IR_Times)])
legend([p1(1) p2(1)],'Kinetic Temperatures','Radiant Temperatures')

%% plotting stream residuals 

stream_res = IRstreamtemps - subibStreamTemps;
mean_r_t = mean(stream_res,1)
mean_r_d = mean(stream_res,2)

%histograms from means
figure; histogram(mean_r_d,6,'FaceColor',[0.5 0.5 0.5]); title('Distribution of Residuals in Space')
xlabel('Temperature Difference (^oC)');ylabel('Frequency')
figure; histogram(mean_r_t,10,'FaceColor',[0.5 0.5 0.5]); title('Distribution of Residuals in Time')
xlabel('Temperature Difference (^oC)');ylabel('Frequency')

%surfaces from mean temps
figure; imagesc(stream_res); ax = gca;
xlabel('Time');ylabel('Control Point')
c = colorbar; c.Label.String ='Radiant Temperature Error (^oC)';colormap('jet')
ax.XTick = [13 49 85 120 156 191 215 251 287 322 357 393 427];
ax.XTickLabel = {'','0:00','','12:00','','0:00','','0:00','','12:00','','0:00',''};

%% F9. plots of residual adjusted temps
load('C:\Users\Emily\Google Drive\Peru_2016\Peru 2016\iButtons_2016\iButton_calibration_correction_16.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\iButton_stream_temps_set1.mat')
load('C:\Users\Emily\Google Drive\IR Codes_Peru 2016\Extracting_control_points\IR_stream_cpts_2016_set1.mat')

stream_res = mean_IRpixeltemps(:,[23:452]) - (subibStreamTemps+corr_16(1:13,:));
cp1 = stream_res(2,:);
cp2 = stream_res(7,:);
cp3 = stream_res(12,:);
cpt_mean_res = (cp1+cp2+cp3)/3;
for i=1:13
adj_IR_T(i,:) = mean_IRpixeltemps(i,23:452) - cpt_mean_res;
end
abs_error = abs(adj_IR_T - (subibStreamTemps+corr_16(1:13,:)));
mean_error = mean(reshape(abs_error,[1,13*430]));

figure; color = colormap(jet)
plot(IR_Times,adj_IR_T(1,:),'Color',[color(1,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(2,:),'Color',[color(6,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(3,:),'Color',[color(11,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(4,:),'Color',[color(17,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(5,:),'Color',[color(22,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(6,:),'Color',[color(27,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(7,:),'Color',[color(32,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(8,:),'Color',[color(37,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(9,:),'Color',[color(43,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(10,:),'Color',[color(48,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(11,:),'Color',[color(54,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(12,:),'Color',[color(59,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(13,:),'Color',[color(64,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,(subibStreamTemps+corr_16(1:13,:)),'k.')
xlabel('Time');ylabel('Temperature (^oC)');
xlim([min(IR_Times) max(IR_Times)])
legend('1-Upstream','2','3','4','5','6','7','8','9','10','11','12','13','Kinetic Temp','Location','northeastoutside')

figure;plot(subibStreamTemps,adj_IR_T,'.')
X=reshape(subibStreamTemps+corr_16(1:13,:),430*13,1);y=reshape(adj_IR_T,430*13,1);
mdl = fitlm(X,y)
%% plots of residual adjusted temps
stream_res = IRstreamtemps - subibStreamTemps;
cp1 = stream_res(4,:);
cp2 = stream_res(8,:);
cp3 = stream_res(13,:);
cpt_mean_res = (cp1+cp2+cp3)/3;
for i=1:13
adj_IR_T(i,:) = IRstreamtemps(i,:) - cpt_mean_res;
end
abs_error = abs(adj_IR_T - subibStreamTemps);
mean_error = mean(reshape(abs_error,[1,13*430]));

figure; color = colormap(jet)
plot(IR_Times,adj_IR_T(1,:),'Color',[color(1,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(2,:),'Color',[color(6,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(3,:),'Color',[color(11,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(4,:),'Color',[color(17,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(5,:),'Color',[color(22,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(6,:),'Color',[color(27,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(7,:),'Color',[color(32,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(8,:),'Color',[color(37,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(9,:),'Color',[color(43,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(10,:),'Color',[color(48,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(11,:),'Color',[color(54,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(12,:),'Color',[color(59,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,adj_IR_T(13,:),'Color',[color(64,:)],'Marker','.','LineStyle','none'); hold on
plot(IR_Times,subibStreamTemps,'k.')
xlabel('Time');ylabel('Temperature (^oC)');
xlim([min(IR_Times) max(IR_Times)])
legend('1-Upstream','2','3','4','5','6','7','8','9','10','11','12','13','Kinetic Temp','Location','northeastoutside')
%% all sets of 3 control points
stream_res = IRstreamtemps - subibStreamTemps;

% Make all combinations of 3 control points
v = 1:1:13; % numbers to choose from
Combos = nchoosek(v,3) % make all cominations

for j = 1:length(Combos)
cp1 = stream_res(Combos(j,1),:);
cp2 = stream_res(Combos(j,2),:);
cp3 = stream_res(Combos(j,3),:);
cpt_mean_res = (cp1+cp2+cp3)/3;

    for i=1:13
    adj_IR_T(i,:) = IRstreamtemps(i,:) - cpt_mean_res;
    end
    
abs_error = abs(adj_IR_T - subibStreamTemps);
mean_error(j) = mean(reshape(abs_error,[1,13*430]));
end

figure; histogram(mean_error,40,'FaceColor',[0.5 0.5 0.5])
xlabel('Mean Absolute Error (^oC)'); ylabel('Frequency')

%% all sets of 5 control points
stream_res = IRstreamtemps - subibStreamTemps;

% Make all combinations of 5 control points
v = 1:1:13; % numbers to choose from
Combos = nchoosek(v,5) % make all cominations

for j = 1:length(Combos)
cp1 = stream_res(Combos(j,1),:);
cp2 = stream_res(Combos(j,2),:);
cp3 = stream_res(Combos(j,3),:);
cp4 = stream_res(Combos(j,4),:);
cp5 = stream_res(Combos(j,5),:);
cpt_mean_res = (cp1+cp2+cp3+cp4+cp5)/5;

    for i=1:13
    adj_IR_T(i,:) = IRstreamtemps(i,:) - cpt_mean_res;
    end
    
abs_error = abs(adj_IR_T - subibStreamTemps);
mean_error(j) = mean(reshape(abs_error,[1,13*430]));
end

figure;histogram(mean_error,40,'FaceColor',[0.5 0.5 0.5])
xlabel('Mean Absolute Error (^oC)'); ylabel('Frequency')