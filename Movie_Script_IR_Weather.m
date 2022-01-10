%% Movie of Data with stream temp and weather data

load('data_reg_Dec1.mat', 'data_reg'); %aligned raw IR image data (3D array)
load('AdjIRtime.mat', 'redoTime'); %time data to plot on x axis
load('IR_Weather_Variables.mat'); %weather variables during IR imaging timeframe
load('raw_IR_temps_stream.mat'); %raw IR stream temps over time at control points
raw_IR_shifted = IRstreamtempstrans+5.5; %shifting IR temps upward to improve fit with iButton sensor temperature data
load('IR_Times.mat')

%the 'true' stream temps that were measured with the ibuttons
ibStreamTemps=xlsread('iButton_stream_temps.xlsx','iButton_stream_temps');
ibStreamTempsTrans=transpose(ibStreamTemps);
subibStreamTemps=ibStreamTemps(:,1:2:end); %every other ibutton data point
%% Movie  
t=char(redoTime); %converting times in datetime format to character formt to put into textbox above graphs to know the time during the movie
v = VideoWriter('IR_Video_Weather'); %creates the file that you want to build the movie in
v.FrameRate = 4; %determines how fast each frame of the movie goes (x frames per second)
v.Quality = 100; %100 to 0 where 100 is highest quality
open(v) %opens the video file so that you can put plots into it

for k = 1:707 %how many images I have; I will plot each image
    f=figure();
    test=data_reg(:,:,k); %takes the first image out of the 3D array
    for i=1:1024; for j=1:768; test2(j,i)=test(j,i,1); end; end; 
  subplot('Position',[0.07 0.5 .53 0.38]) %plot of IR picture
    image(test2,'CDataMapping','scaled')
    colormap(jet) %sets color scheme for the temp
      caxis([-10 30]) %set range of z axis (temp)
      xmin = min(0); %set range of x axis
      xmax = max(1024);
      xlim([xmin xmax]);
      ymin = min(0); %set range of y axis
      ymax = max(768);
      ylim([ymin ymax]);
  str=t(k,:); %assigns the date and time of each image into a textbox above the image plot
  dim=[.32 .69 .3 .3]; %sets dimensions of textbox
  annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','non','FontWeight','bold','FontSize',15); %puts the date and time into the textbox and sets its properties
      c=colorbar; %puts in the legend for the temperature gradient (z axis)
      c.Label.String = 'Pixel Temperature (^oC)'
      c.Label.FontWeight ='bold';
      
   subplot('Position',[0.08 0.1 .5 .3]) %stream temperature verses time
      plot(IR_Times, raw_IR_shifted,IR_Times,subibStreamTemps); %plot IR temps and ibutton temps
      hold on
      plot(IR_Times(k,:),raw_IR_shifted(k,:),'r.') %plot just the specific point in time that you are tracking in the image
      datetick('x', 'mm/dd HHAM','keepticks'); 
    set (gca,'FontSize', 9);
    xmin = min(736165.65);
    xmax = max(736170.6);
    xlim([xmin xmax]);
    ymin = min(-5);
    ymax = max(30);
    ylim([ymin ymax]);
    xlabel('Time','FontWeight','bold')
    ylabel('Temperature (^oC)','FontWeight','bold')  
    
  subplot('Position',[0.75 0.71 .22 .21]) %temp versus time
      plot(IR_Times, Station_Temp); %plot air temp
      hold on
      plot(IR_Times(k,:),Station_Temp(k,:),'r.') %plot just the specific point in time that you are tracking in the image
      datetick('x', 'mm/dd HHAM','keepticks'); 
    set (gca,'FontSize', 9);
    xmin = min(736165.65);
    xmax = max(736170.6);
    xlim([xmin xmax]);
    ymin = min(-5);
    ymax = max(20);
    ylim([ymin ymax]);
    ylabel('Air Temp (^oC)','FontWeight','bold')  
    
    subplot('Position',[0.75 0.48 .22 .21]) %RH versus time
      plot(IR_Times, Station_RH); %plot RH
      hold on
      plot(IR_Times(k,:),Station_RH(k,:),'r.') %plot just the specific point in time that you are tracking in the image
      datetick('x', 'mm/dd HHAM','keepticks'); 
    set (gca,'FontSize', 9);
    xmin = min(736165.65);
    xmax = max(736170.6);
    xlim([xmin xmax]);
    ymin = min(10);
    ymax = max(100);
    ylim([ymin ymax]);
    ylabel('RH (%)','FontWeight','bold')  
    
    subplot('Position',[0.75 0.252 .22 .21]) %Solar versus time
      plot(IR_Times, Solar); %plot solar
      hold on
      plot(IR_Times(k,:),Solar(k,:),'r.') %plot just the specific point in time that you are tracking in the image
      datetick('x', 'mm/dd HHAM','keepticks'); 
    set (gca,'FontSize', 9);
    xmin = min(736165.65);
    xmax = max(736170.6);
    xlim([xmin xmax]);
    ymin = min(0);
    ymax = max(1300);
    ylim([ymin ymax]);
    ylabel('Solar (W/m^2)','FontWeight','bold')  
    
      
    subplot('Position',[0.75 0.02 .22 .21]) %wind versus time
      plot(IR_Times, Wind); %plot wind
      hold on
      plot(IR_Times(k,:),Wind(k,:),'r.') %plot just the specific point in time that you are tracking in the image
      datetick('x', 'mm/dd HHAM','keepticks'); 
    set (gca,'FontSize', 9);
    xmin = min(736165.65);
    xmax = max(736170.6);
    xlim([xmin xmax]);
    ymin = min(0);
    ymax = max(20);
    ylim([ymin ymax]);
    xlabel('Time','FontWeight','bold')
    ylabel('Wind (m/s)','FontWeight','bold')  
   
    drawnow 
      pause(0.1)
      F=getframe(f); %captures the figure
  writeVideo(v,F);
  close all %closes the figure so that I don't run out of memory
end  
close(v);
