%This code runs the filtering and smoothing on the raw strain data

angles=[36.00000	271.00000	178.00000	80.00000	280.00000	182.00000	202.00000	274.0000 ];
m=     [-12.19000	-11.07000	-10.07000	-10.20000	-11.08000	-14.06000	-10.24000	-9.65000 ];
arms=  [130.00000	130.00000	130.00000	130.00000	130.00000	130.00000	130.00000	130.0000 ];

%Import_SG_Info
%angles=SG_Info.Angle;
%m=SG_Info.m;
%arms=SG_Info.armsmm;


%% 1. Filter and convert from raw data to strain 
load('T18_21_AllData.mat')
data=T18_21_AllData;  
%plot(datetime(T18_21_AllData(5000000:end-100,1), 'ConvertFrom', 'datenum'),T18_21_AllData(5000000:end-100,9));

window=5000;
for col=2:9
    col_data=data(:,col);

    %Filtering for outliers and NaN's
    col_data(find(col_data>=100))=0; %positive outliers
    col_data(find(col_data<=-100))=0; %negative outliers
    col_data(find(isnan(col_data==1)))=0; %NaN's
    
    %Smoothing with a 500 point running mode
    col_data=Running_mode(col_data,window);

    %Calibrate with arm length and voltage to extension gradient 'm'
    col_data=col_data/(arms(col-1)*m(col-1));
    if col==2
        data_out=col_data;
    else
        data_out=cat(2,data_out,col_data);% cat together the infor from all the pairs
    end
end % end loop over cols
data=cat(2,data(:,1),data_out);
T_cal=data; 
%save('T1_cal','T1_cal')   
%%
plot(datetime(T_cal(6000000:end-100,1), 'ConvertFrom', 'datenum'),T_cal(6000000:end-100,9));

%% zoom in plots
sample=2e5:2.2e6; %
%sample=50000:size(T18_21_AllData,1)-100;
%sample=6000000:6010000;
%sample=6008000:6010000;

plot(datetime(data(sample,1), 'ConvertFrom', 'datenum'),data(sample,9));
hold on
plot(datetime(data(sample,1), 'ConvertFrom', 'datenum'),Running_mode1(data(sample,9),500));
plot(datetime(data(sample,1), 'ConvertFrom', 'datenum'),Running_mode1(data(sample,9),5000));
plot(datetime(data(sample,1), 'ConvertFrom', 'datenum'),Running_mode1(data(sample,9),10000));
legend('Raw data', '1000 point (4 min) mode','10000 point (40 min) mode','20000 point (80 min) mode')
%%


load('C:\Users\Toby\Documents\MATLAB\Wytham_Summary\Wind_data\all_my_wind1.mat')
rows_T1=find(data(:,1)>data(1.5e5,1) & data(:,1)<data(2e6,1));
rows_wind=find(all_my_wind1(:,1)>data(1.5e5,1) & all_my_wind1(:,1)<data(2e6,1));
%subplot(2,1,1)
plot(datetime(data(rows_T1,1), 'ConvertFrom', 'datenum'),data(rows_T1,[2])*1000);
hold on
plot(datetime(all_my_wind1(rows_wind,1), 'ConvertFrom', 'datenum'),all_my_wind1(rows_wind,3));

%% Overlay simple and robust maxima finders
%This part of the data is continuous and error free, so I can just loop by index.
%I plot each hour and select the maxima using both methds
hold on
temp_data1=data(rows_T1(1:691200+1),[2]);
temp_time=datetime(data(rows_T1(1:691200+1),1), 'ConvertFrom', 'datenum');
colors=gray(10);
h1=plot(temp_time,temp_data1*1000,'Color',colors(6,:));
hold on
colors=winter(10);
h2=plot(datetime(all_my_wind1(rows_wind(1:17280+1),1), 'ConvertFrom', 'datenum'),all_my_wind1(rows_wind(1:17280+1),3),'Color',colors(6,:));
for i=1:48 %hourly samples
    start_row=(i-1)*14400+1; end_row=i*14400;
    
    [M I]=max(temp_data1(start_row:end_row));
    %There could be multiple of the same max - so I just need to plot
    plot(temp_time(I+start_row),M*ones(length(I),1)*1000,'rx')
    %straight off - but then what about the datenum
    %row_number_simple=start_row+find
    
    temp_robust=robust_max(temp_data1(start_row:end_row));
    plot([temp_time(start_row) temp_time(end_row)],1000*[temp_robust temp_robust],'Color','green')
    plot([temp_time(end_row) temp_time(end_row)],[0 7],'Color','black')
end

legend('Strain data', 'Wind data', 'Simple max', 'Robust max', 'Location', 'NorthWest')
legend boxoff
title('Hourly maximum strain selection')

%% 2. Convert T1_cal strain data to max strain and NE
%load('T1_cal.mat')
data=T_cal;  
data_out=cat(3,data(:,1),data(:,1));

tot_cols=0;
this_col=0;
for pair=2:5

    this_col=this_col+2;
    first=(pair-1)*2; second=first+1;
    col_data=data(:,[first second]);

    for col=1:2 % DO I REALLY NEED TO FILER AGAIN?
        %Filtering for outliers and NaN's
        col_data(find(col_data(:,col)>=1e-2),col)=0; %positive outliers
        col_data(find(col_data(:,col)<=-1e-2),col)=0; %negative outliers
        col_data(find(isnan(col_data(:,col)==1)),col)=0; %NaN's
    end

    [max_strain, angle] = Solve_vector( col_data(:,1), col_data(:,2), angles(this_col-1), angles(this_col) );
    data_out=cat(2,data_out,cat(3,max_strain,angle));% cat together the infor from all the pairs
end % end loop over pairs
T_MaxStrain=data_out; 
%save('T1_MaxStrain','T1_MaxStrain')
     
% NE section
%load('T1_MaxStrain.mat')
data=T_MaxStrain;  

for col=1:5
    N_strain(:,col)=data(:,col,1).*cosd(data(:,col,2));
    E_strain(:,col)=data(:,col,1).*sind(data(:,col,2));
end
data_out=cat(3,N_strain,E_strain);
data_out(:,1,:)=data(:,1,:);
T_NE=data_out; 
%save('T1_NE','T1_NE')
    

%% plot samples
sample=6008000:6010000;
subplot(3,1,1)
plot(datetime(T_NE(sample,1,1), 'ConvertFrom', 'datenum'),T_NE(sample,5,1));


title('Strain on North facing side of trunk')
ylabel('strain')
subplot(3,1,2)
plot(datetime(T_NE(sample,1,1), 'ConvertFrom', 'datenum'),T_NE(sample,5,2));

title('Strain on East facing side of trunk')
ylabel('strain')
subplot(3,1,3)
plot(datetime(T_MaxStrain(sample,1,1), 'ConvertFrom', 'datenum'),T_MaxStrain(sample,5,1));
title('Max Strain format')
ylabel('strain')
pause
close all

%% 3. Loop over hourly wind data and select max strain in each hour

load('hourly_ECN_Walkway_wind.mat')
wind=hourly_ECN_Walkway_wind;
clearvars hourly_ECN_Walkway_wind
data_in=T_MaxStrain(:,:,1);   

% Averaging process starts here
no_cols=size(data_in,2);  
data_out_simple_max=NaN(5951,size(data_in,2));
data_out_robust_max=NaN(5951,size(data_in,2));

%Now we loop over the hours of ECN data and find the max of tree strain data within that hour
for hour= 2:5951 
    hour
    low=wind(hour-1,1); high=wind(hour,1); %Select the datenums at the start and end of the hour
    rows=find(data_in(:,1)>=low & data_in(:,1)<high); %Find strain data within that hour
    if length(rows)<14000 continue %If there is too little strain data skip it
    end
    if wind(hour,3)==0 continue  %If there was no wind (or more often errors gives 0) skip it
    end
    length(rows)
    data_out_simple_max(hour,:)=max(data_in(rows,:)); %simple data out - this should be good enough.
    data_out_robust_max(hour,:)= robust_max(data_in(rows,:) );
    %plot(data_in(rows,2:end))
 
end

T_hourly_simple_max=cat(2,data_out_simple_max(:,1),wind(:,2:3),data_out_simple_max(:,2:end));
T_hourly_robust_max=cat(2,data_out_robust_max(:,1),wind(:,2:3),data_out_robust_max(:,2:end));
scatter(T_hourly_simple_max(:,2),T_hourly_simple_max(:,5))
hold on
scatter(T_hourly_robust_max(:,2),T_hourly_robust_max(:,5))
ylabel('Strain')
xlabel('Wind')
legend('Simple maxima' , 'Modal maxima','Location','NorthWest')
legend boxoff
pause


%% 4. Calculate critical wind speed - hourly

breaking_strain=6.86e-3; %this is taken from the literature and up for debate.
CWS_power_winter = zeros(4,3);
FITS_power_winter = zeros(4,4);
CWS_linear_winter = zeros(4,3);
FITS_linear_winter = zeros(4,3);

%Loops over each tree 
for tree=1:4 
    tree
    col=tree+3;
    %Imports ECN hourly data and strain data %2=ECN mean 7=canopy max
    Winter = cat(2,T_hourly_simple_max(1:2999,2),T_hourly_simple_max(1:2999,col)) ;
    Winter(Winter(:,2)==-inf,2)=NaN;     %Checks for infinities
    for col = 1:2
        %Removes the NaN
        Winter = Winter(isnan(Winter(:,col))== 0,:);
    end
  
    strain=Winter(:,2);
    wind=Winter(:,1);
    % Winter
   %Function: fitpower bisquare. Green Properties pulls the strain
   %specific to the tree from the excel            %strain      %wind
   [fit_winter1, gof_winter1] = fit_power_bisquare(strain, wind,breaking_strain);
   CI = confint(fit_winter1);        %Calculates the upper and lower error bounds using the confint function
   CWS_power_winter(tree,1:3) = [(fit_winter1.a)*(breaking_strain^fit_winter1.b) (CI(1,1))*(breaking_strain^CI(1,2)) (CI(2,1))*(breaking_strain^CI(2,2))];
   %Makes tables for R^2, root mean square error, and function parameters(A,B)
   FITS_power_winter(tree,1:4) = [(gof_winter1.rsquare) (gof_winter1.rmse) (fit_winter1.a) (fit_winter1.b)];

   %Function: fitpower force to 0. Green Properties pulls the strain
   %specific to the tree from the excel
   [fit_winter2, gof_winter2] = fit_bisquare_FORCE0(strain,wind);
   CI = confint(fit_winter2);%Calculates the upper and lower error bounds using the confint function
   CWS_linear_winter(tree,1:3) = [(fit_winter2.p1)*sqrt(breaking_strain) (CI(1,1))*sqrt(breaking_strain) (CI(2,1))*sqrt(breaking_strain)];
   %Makes tables for R^2, root mean square error, and function
   %parameters(A,B)         
   FITS_linear_winter(tree,1:3) = [(gof_winter2.rsquare) (gof_winter2.rmse) (fit_winter2.p1)];


    % Plot of Power bisquare

  
    subplot(1,2,1)
    color=jet(10);
    h1 = plot( strain, wind,'+','Color',color(4,:));
    h1fit=line(linspace(0,breaking_strain,100),fit_winter1.a*(linspace(0,breaking_strain,100).^(fit_winter1.b)));
    %h1fit=line(linspace(0,breaking_strain.^(fit_winter1.b),100).^(fit_winter1.b),...
    %    fit_winter1.a*(linspace(0,breaking_strain.^(fit_winter1.b),100).^(fit_winter1.b)));
    xlim([0 breaking_strain])  %Set axis for extrapolation
    hold on
    h2fit=plot(fit_winter1,'predobs');
    plot(breaking_strain*ones(1000,1),linspace(0,40,1000)','r--')
    set([h1fit],'Color',color(4,:))
    set([h2fit],'Color',color(4,:))
    legend off
    legend( [h1],'Winter' , 'Location', 'NorthWest' );
    legend boxoff

    grid off
    xlim([0 breaking_strain+0.002])
    ylim([0 CWS_power_winter(tree,1)+5])
    title('CWS - variable exponent')
    ylabel('ECN Mean Wind Speed')
    xlabel('nth root of strain')
    exponent=num2str(round(fit_winter1.b,2));
    %textLoc({['Winter CWS = ' num2str(CWS_power_winter(tree,1))],['Wind = ' num2str(round(fit_winter1.a)) ' strain ^' exponent(1) '^' exponent(2) '^' exponent(3) '^' exponent(4)] },'SouthEast');
    %str=strcat('Fit_Images/CWS Mean ECN Tree  ', num2str(tree));

    subplot(1,2,2)
    winter_residuals1=wind-fit_winter1.a.*(strain.^(fit_winter1.b));    
    histogram(winter_residuals1,'facecolor',color(4,:),'facealpha',.5,'edgecolor',color(4,:));
    title residuals
    pause
    close all
    
    FORCED_FIT=0;
    if FORCED_FIT==1
        % Plot of Forced fit
        subplot(1,2,1)
        h1 = plot( strain, wind,'+','Color',color(4,:));
        xlim([0 breaking_strain])  %Set axis for extrapolation
        hold on
        plot(breaking_strain*ones(1000,1),linspace(0,50,1000)','r--')
        %h1fit=plot(fit_winter2,'predobs');
        h1fit=line(linspace(0,breaking_strain,100),fit_winter2.p1*(linspace(0,breaking_strain,100)).^(1/2));

        set(h1fit,'Color',color(4,:))
        legend off
        legend( [h1],'Winter' , 'Location', 'NorthWest' );
        legend boxoff
        grid off
         xlim([0 breaking_strain+0.002])
        ylim([0 CWS_linear_winter(tree,1)+5])
        xlabel sqrt(Strain)
        ylabel('ECN Mean Wind Speed')
        title('CWS - exponent == 2')
        str=strcat('Fit_Images/CWS Mean ECN Tree  ', num2str(tree));
        textLoc({['Winter CWS = ' num2str(CWS_linear_winter(tree,1))],['Wind = ' num2str(round(fit_winter2.p1)) 'strain ^0^.^5']},'SouthEast');
        
        subplot(1,2,2)
        winter_residuals2=wind-fit_winter2.p1*(strain.^(0.5));     
        histogram(winter_residuals2,'facecolor',color(4,:),'facealpha',.5,'edgecolor',color(4,:));
        title residuals
        pause
    end


    close all
end

CWS=cat(2,CWS_power_winter, CWS_linear_winter);
FITS=cat(2,FITS_power_winter, FITS_linear_winter);
mat2clip([CWS(1, [1 4]) CWS(2, [1 4]) CWS(3, [1 4]) CWS(4, [1 4])])


mat2clip([ CWS(1:4, 4)'; CWS(1:4, 1)'])


mat2clip([CWS(1:4, 1)'])





