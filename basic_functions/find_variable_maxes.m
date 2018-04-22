function [ output ] = find_variable_maxes( MaxStrain, wind_data, time_average, wind_points_threshold,strain_points_threshold )

%% Get rid of infinities
for col=1:size(MaxStrain,2)
    MaxStrain(isinf(MaxStrain(:,col)==1),col)=NaN;
end
for col=1:size(wind_data,2)
    wind_data(isinf(wind_data(:,col)==1),col)=NaN;
end

%% Split the data by the datestamp
[ splits ] = find_splits( MaxStrain(:,:,1), 0.006945 ); %Find the splits in the strain data 
count_too_little=0; count_too_much=0; count_sections=0;

for section=1:length(splits)-1 %Loop over these splits
    strain_section=MaxStrain(splits(section):splits(section+1)-1,:,1); %Define temp strain vector between two splits
    [ wind_data_out ] = find_wind_overlap_continuous( strain_section, wind_data ); %Find concurrent wind data
    if isnan(wind_data_out)==1; continue
    end
        
    len=(wind_data_out(end,1)-wind_data_out(1,1))/time_average
    %linspace floor
        temp_mean_max=nan(floor(len),size(wind_data_out,2)+size(strain_section,2)-1);
    for i=1:floor(len) %this is the loop over hours / seconds / minutes
        select_time(i)=wind_data_out(1,1)+(i)*time_average; %Hour ending format
        select_wind=wind_data_out(find(wind_data_out(:,1)>=wind_data_out(1,1)+(i-1)*time_average & wind_data_out(:,1)<=wind_data_out(1,1)+i*time_average),[2:size(wind_data_out,2)]);
        if length(select_wind)<wind_points_threshold; continue
        end
        select_strain=strain_section(find(strain_section(:,1)>=wind_data_out(1,1)+(i-1)*time_average & strain_section(:,1)<=wind_data_out(1,1)+i*time_average),[2:size(strain_section,2)]); %Find 10s strain data
        if length(select_strain)<strain_points_threshold; continue
        end
        temp_mean_max(i,[1:size(wind_data_out,2)+size(strain_section,2)-1])=cat(2,select_time(i),mean(select_wind), max(select_strain));
        
    end
    count_sections=count_sections+1;
    [count_sections length(splits)]
    if count_sections==1
        output=temp_mean_max;
    else
        output=cat(1,output,temp_mean_max);
    end
    clearvars temp_mean_max
end


end %end of function

