%This code runs the filtering and smoothing on the raw strain data
%The calibration information for T18_21 are:

angles=[220.00000	117.00000	200.00000	310.00000	288.00000	190.00000	190.00000	278.00000];
m= [-12.04156	-13.17500	-14.46494	-15.67004	-12.65817	-18.40704	-12.59363	-10.38261];
arms=  [130.00000	130.00000	130.00000	130.00000	130.00000	130.00000	130.00000	130.0000 ];


%% 1. Filter and convert from raw data to strain 
load('C:\Users\Toby\Dropbox\Tobys_Stuff\MATLAB\Wytham_Summary\Strain_data\AllData_Files\T18_21_AllData.mat')
data=T18_21_AllData;  
%dlmwrite('T18_21_AllData.txt',T18_21_AllData(:,:,1),'precision',13)

for col=2:9
    col_data=data(:,col);

    %Filtering for outliers and NaN's
    col_data(find(col_data>=100))=0; %positive outliers
    col_data(find(col_data<=-100))=0; %negative outliers
    col_data(find(isnan(col_data==1)))=0; %NaN's
    
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
  
%% Calculate NE (uploaded version) 

for col=2:2:8
    unsmoothed_NE(:,col)=  T_cal(:,col).*cosd(angles(col-1))+T_cal(:,col+1).*sind(angles(col));
    unsmoothed_NE(:,col+1)=T_cal(:,col).*sind(angles(col-1))+T_cal(:,col+1).*cosd(angles(col));
end
unsmoothed_NE(:,1)=T_cal(:,1);
%T18_21_NE_unsmoothed=unsmoothed_NE;
%dlmwrite('T18_21_NE_unsmoothed.txt',T18_21_NE_unsmoothed,'precision',13)

%% ==============================================================
%%          START FROM HERE IF USING NE UNSMOOTHED DATA
%% ==============================================================

%% Smooth using a running mode
window_size=5000;
for col=2:9
    NE_smoothed(:,col)=Running_mode(T18_21_NE_unsmoothed(:,col),window_size);
end
NE_smoothed(:,1)=T18_21_NE_unsmoothed(:,1);

%% 2. Convert T1_cal strain data to max strain and NE
data=NE_smoothed;  

T_MaxStrain(:,1)=data(:,1); 
T_MaxStrain(:,2)=sqrt(data(:,2).^2+data(:,3).^2); 
T_MaxStrain(:,3)=sqrt(data(:,4).^2+data(:,5).^2); 
T_MaxStrain(:,4)=sqrt(data(:,5).^2+data(:,7).^2); 
T_MaxStrain(:,5)=sqrt(data(:,6).^2+data(:,9).^2); 

%% Plot two channels and max strain

sample=6008000:6010000;
subplot(3,1,1)
plot(datetime(data(sample,1), 'ConvertFrom', 'datenum'),data(sample,8));
title('Strain on North facing side of trunk')
ylabel('strain')

subplot(3,1,2)
plot(datetime(data(sample,1), 'ConvertFrom', 'datenum'),data(sample,9));
title('Strain on East facing side of trunk')
ylabel('strain')

subplot(3,1,3)
plot(datetime(T_MaxStrain(sample,1), 'ConvertFrom', 'datenum'),T_MaxStrain(sample,5));
title('Max Strain')
ylabel('strain')


%% 3. Loop over hourly wind data and select max strain in each hour

load('hourly_ECN_Walkway_wind.mat')
wind=hourly_ECN_Walkway_wind;
clearvars hourly_ECN_Walkway_wind
data_in=T_MaxStrain;   

% Averaging process starts here
no_cols=size(data_in,2);  
data_out=NaN(5951,size(data_in,2));

%Now we loop over the hours of ECN data and find the max of tree strain data within that hour
for hour= 2:5951 
    hour
    low=wind(hour-1,1); high=wind(hour,1); %Select the datenums at the start and end of the hour
    rows=find(data_in(:,1)>=low & data_in(:,1)<high); %Find strain data within that hour
    if length(rows)<14000 continue %If there is too little strain data skip it
    end
    if wind(hour,3)==0 continue  %If there was no wind (or more often errors give 0) skip it
    end
    length(rows)
    data_out(hour,:)=max(data_in(rows,:)); %simple data out - this should be good enough.
end

T_hourly=cat(2,data_out(:,1),wind(:,2:3),data_out(:,2:end));


%% 4. Calculate critical wind speed - hourly

breaking_strain=5e-3; %this is taken from the literature and up for debate.
%Loops over each tree 
for tree=1:4 
    tree
    subplot(2,2,tree)
    col=tree+3;
    %Imports ECN hourly data and strain data %2=ECN mean 3=ECN max
    Winter = cat(2,T_hourly(1:2999,3),T_hourly(1:2999,col)) ;
    Winter(Winter(:,2)==-inf,2)=NaN;     %Checks for infinities
    for col = 1:2
        %Removes the NaN
        Winter = Winter(isnan(Winter(:,col))== 0,:);
    end
    strain=Winter(:,2);
    wind=Winter(:,1);

    color=jet(10);
    ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    [fit_square_data, gof_square_data] = fit(wind.^2,strain,ft,opts);
    save_norms(tree)=fit_square_data.a;
    CWS_data(tree) = [(breaking_strain/(save_norms(tree))).^(1/2)]; 
       
    h1 = plot( wind.^2, strain,'+','Color',color(2,:));
    hold on
    h_data_fit=line(linspace(0,25,100).^2,(save_norms(tree))*(linspace(0,25,100).^2),'Color',color(5,:));
    hbreak=plot(linspace(0,CWS_data(tree).^2,1000)',breaking_strain*ones(1000,1),'r--');
    ylim([1e-5 breaking_strain+0.002])
    xlabel('Max ECN wind speed')
    ylabel('Max hourly strain')
    pause
  
end






