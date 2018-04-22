function [ wind_data_out ] = find_wind_overlap_continuous( strain_data_in, wind_data_in )
   %Find wind sample overlapping with my strain sample, assuming sample is continuous
   
    time=strain_data_in(isnan(strain_data_in(:,1))==0,1);
    if length(time)>0;
 
        find_wind_overlap=find(wind_data_in(:,1)>=time(1,1) & wind_data_in(:,1)<=time(length(time),1));
        wind_data_out=wind_data_in(find_wind_overlap,:);

        if length(wind_data_out)==0
            wind_data_out=NaN;
            disp('NO WIND DATA HERE')
        end
    else
        wind_data_out=NaN;
        disp('NO NON NAN DATA')
    end   
    
end

