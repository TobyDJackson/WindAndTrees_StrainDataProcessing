function [ where_splits ] = find_splits( data_in, gap )


% Set up empty vectors
diff=nan(length(data_in),1);
splits=zeros(length(data_in),1); splits(1)=1; splits(end)=1;
hour=0.04167; tenmin=0.006945; 

%Calculate difference between consectutive datenums, and find splits
for i=2:length(data_in)
    diff(i)=data_in(i,1)-data_in(i-1,1);
    if diff(i)>=gap %this represents about a 30 minute gap
        splits(i)=1;
    end
end
where_splits=find(splits==1);

%Loop over split start-end point pairs - to get an overview
%row=0;
%for j=2:2:length(where_splits)
%    row=row+1; %this serves as as counter since we are looping in 2's
%    split_lengths(row)=(where_splits(j)-where_splits(j-1))/(4*60*60); %This should give the split length in hours, assuming 4Hz continuous
%    check_lengths(row)=(data_in(where_splits(j)-1,1)-data_in(where_splits(j-1),1))/hour ; %Same result, based on datenum
%    first_part(row)=datetime(data_in(where_splits(j-1),1), 'ConvertFrom', 'datenum'); %Save start datetimes
%end

end %end of fn

