%This code runs the filtering and smoothing on the raw strain data

angles=[36.00000	271.00000	178.00000	80.00000	280.00000	182.00000	202.00000	274.0000 ];
m=     [-12.19000	-11.07000	-10.07000	-10.20000	-11.08000	-14.06000	-10.24000	-9.65000 ];
arms=  [130.00000	130.00000	130.00000	130.00000	130.00000	130.00000	130.00000	130.0000 ];

%Import_SG_Info
%angles=SG_Info.Angle;
%m=SG_Info.m;
%arms=SG_Info.armsmm;


%% 1. Filter and convert from raw data to strain 
load('T1_AllData.mat')
data=T1_AllData;  
%plot(datetime(T1_AllData(5000000:end-100,1), 'ConvertFrom', 'datenum'),T1_AllData(5000000:end-100,9));

window=50000;
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
T1_cal=data; 
%save('T1_cal','T1_cal')   
%plot(datetime(T1_cal(5000000:end-100,1), 'ConvertFrom', 'datenum'),T1_cal(5000000:end-100,9));

%% zoom in plots
%sample=5000000:size(T1_AllData,1)-100;
%sample=6000000:6010000;
sample=6008000:6010000;
subplot(2,1,1)
plot(datetime(T1_AllData(sample,1), 'ConvertFrom', 'datenum'),T1_AllData(sample,9));
title('raw data')
ylabel('mV')
subplot(2,1,2)
plot(datetime(T1_cal(sample,1), 'ConvertFrom', 'datenum'),T1_cal(sample,9));
title('calibrated data')
ylabel('strain')
pause
close all


%% 2. Convert T1_cal strain data to max strain and NE
%load('T1_cal.mat')
data=T1_cal;  
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
T1_MaxStrain=data_out; 
%save('T1_MaxStrain','T1_MaxStrain')
     
% NE section
%load('T1_MaxStrain.mat')
data=T1_MaxStrain;  

for col=1:5
    N_strain(:,col)=data(:,col,1).*cosd(data(:,col,2));
    E_strain(:,col)=data(:,col,1).*sind(data(:,col,2));
end
data_out=cat(3,N_strain,E_strain);
data_out(:,1,:)=data(:,1,:);
T1_NE=data_out; 
%save('T1_NE','T1_NE')
    

%% plot samples
sample=6008000:6010000;
subplot(3,1,1)
plot(datetime(T1_NE(sample,1,1), 'ConvertFrom', 'datenum'),T1_NE(sample,5,1));
title('Strain on North facing side of trunk')
ylabel('strain')
subplot(3,1,2)
plot(datetime(T1_NE(sample,1,1), 'ConvertFrom', 'datenum'),T1_NE(sample,5,2));
title('Strain on East facing side of trunk')
ylabel('strain')
subplot(3,1,3)
plot(datetime(T1_MaxStrain(sample,1,1), 'ConvertFrom', 'datenum'),T1_MaxStrain(sample,5,1));
title('Max Strain format')
ylabel('strain')
pause
close all

%% 3. Loop over hourly wind data and select max strain in each hour

load('hourly_ECN_Walkway_wind.mat')
wind=hourly_ECN_Walkway_wind;
clearvars hourly_ECN_Walkway_wind
data_in=T1_MaxStrain(:,:,1);   

% Averaging process starts here
no_cols=size(data_in,2);  
data_out_simple_max=NaN(5951,size(data_in,2));

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
    plot(data_in(rows,2:end))
 
end

T1_hourly_simple_max=cat(2,data_out_simple_max(:,1),wind(:,2:3),data_out_simple_max(:,2:end));
scatter(T1_hourly_simple_max(:,2),T1_hourly_simple_max(:,5))
xlabel('Strain')
ylabel('Wind')




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
    Winter = cat(2,T1_hourly_simple_max(1:2999,2),T1_hourly_simple_max(1:2999,col)) ;
    Winter(Winter(:,2)==-inf,2)=NaN;     %Checks for infinities
    for col = 1:2
        %Removes the NaN
        Winter = Winter(isnan(Winter(:,col))== 0,:);
    end
  
    x=Winter(:,2);
    y=Winter(:,1);
    % Winter
   %Function: fitpower bisquare. Green Properties pulls the strain
   %specific to the tree from the excel            %strain      %wind
   [fit_winter1, gof_winter1] = fit_power_bisquare(Winter(:,2), Winter(:,1),breaking_strain);
   CI = confint(fit_winter1);        %Calculates the upper and lower error bounds using the confint function
   CWS_power_winter(tree,1:3) = [(fit_winter1.a)*(breaking_strain^fit_winter1.b) (CI(1,1))*(breaking_strain^CI(1,2)) (CI(2,1))*(breaking_strain^CI(2,2))];
   %Makes tables for R^2, root mean square error, and function parameters(A,B)
   FITS_power_winter(tree,1:4) = [(gof_winter1.rsquare) (gof_winter1.rmse) (fit_winter1.a) (fit_winter1.b)];

   %Function: fitpower force to 0. Green Properties pulls the strain
   %specific to the tree from the excel
   [fit_winter2, gof_winter2] = fit_bisquare_FORCE0(Winter(:,2), Winter(:,1), sqrt(breaking_strain));
   CI = confint(fit_winter2);%Calculates the upper and lower error bounds using the confint function
   CWS_linear_winter(tree,1:3) = [(fit_winter2.p1)*sqrt(breaking_strain) (CI(1,1))*sqrt(breaking_strain) (CI(2,1))*sqrt(breaking_strain)];
   %Makes tables for R^2, root mean square error, and function
   %parameters(A,B)         
   FITS_linear_winter(tree,1:3) = [(gof_winter2.rsquare) (gof_winter2.rmse) (fit_winter2.p1)];


    % Plot of Power bisquare

  
    subplot(2,2,1)
    color=jet(10);
    h1 = plot( (Winter(:,2).^(fit_winter1.b)), Winter(:,1),'+','Color',color(4,:));
    h1fit=line(linspace(0,breaking_strain.^(fit_winter1.b),100).^(fit_winter1.b),fit_winter1.a*(linspace(0,breaking_strain.^(fit_winter1.b),100).^(fit_winter1.b)));
    xlim([0 (breaking_strain.^(fit_winter1.b))])  %Set axis for extrapolation
    hold on
    %h1fit=line(fit_winter2,'predobs');
    plot((breaking_strain.^(fit_winter1.b))*ones(1000,1),linspace(0,25,1000)','r--')
    set(h1fit,'Color',color(4,:))
    legend off
    legend( [h1],'Winter' , 'Location', 'NorthWest' );
    legend boxoff

    grid off
    xlim([0 breaking_strain.^(fit_winter1.b)+0.02])
    ylim([0 CWS_power_winter(tree,1)+5])
    title('CWS - variable exponent')
    ylabel('ECN Mean Wind Speed')
    xlabel('nth root of strain')
    textLoc(['Winter CWS = ' num2str(CWS_power_winter(tree,1))],'SouthEast');
    %str=strcat('Fit_Images/CWS Mean ECN Tree  ', num2str(tree));

    subplot(2,2,3)
    winter_residuals1=Winter(:,1)-fit_winter1.a.*(Winter(:,2).^(fit_winter1.b));    
    histogram(winter_residuals1,'facecolor',color(4,:),'facealpha',.5,'edgecolor',color(4,:));
    title residuals
    
    % Plot of Forced fit
    subplot(2,2,2)
    h1 = plot( sqrt(Winter(:,2)), Winter(:,1),'+','Color',color(4,:));
    xlim([0 sqrt(breaking_strain)])  %Set axis for extrapolation
    hold on
    h1fit=plot(fit_winter2,'predobs');
    plot(sqrt(breaking_strain)*ones(1000,1),linspace(0,25,1000)','r--')
    set(h1fit,'Color',color(4,:))
    legend off
    legend( [h1],'Winter' , 'Location', 'NorthWest' );
    legend boxoff
    grid off
    xlim([0 0.1])
    ylim([0 inf])
    xlabel sqrt(Strain)
    ylabel('ECN Mean Wind Speed')
    title('CWS - exponent == 2')
    str=strcat('Fit_Images/CWS Mean ECN Tree  ', num2str(tree));
    textLoc(['Winter CWS = ' num2str(CWS_linear_winter(tree,1))],'SouthEast');


    subplot(2,2,4)
    winter_residuals2=Winter(:,1)-fit_winter2.p1*(Winter(:,2).^(0.5));     
    histogram(winter_residuals2,'facecolor',color(4,:),'facealpha',.5,'edgecolor',color(4,:));
    title residuals
    pause
    close all
end

CWS=cat(2,CWS_power_winter, CWS_linear_winter);
FITS=cat(2,FITS_power_winter, FITS_linear_winter);
mat2clip(CWS(:, [1 4]))










