function [ modal_maxima ] = robust_max( hour_data )
% Re-creating the robust max from Barry Gardiners description
% Should work on all columns of a matrix.
%%
%1.       Break the period of interest (normally 60 minutes) into 20 periods (normally each of 3 minutes).
%2.       Find the absolute maximum strain (or wind speed or bending moment) within each separate period. This gives you 20 individual values.
temp_data=nan(20,size(hour_data,2));
len=size(hour_data,1);
for i=1:20
    temp_data(i,:)=max(hour_data((i-1)*len/20+1:i*len/20,:));
end
%%
%3.       Rank the values (m) from 1 to 20 with 1 being the minimum and 20 being the maximum.
sort_data=sort(temp_data);
%%
%4.       Calculate the reduce variate of the CDF: y = -ln(-ln(m/(N+1)) 
%where m is the ranking (1, 2, 3, .., 20) and N = 20 (number of samples).
a=1:20;
reduced_variate=-1*log(-1*log(a/21))';
%%
%5.       Plot your values (strain, wind speed, bending moment, etc.) against the reduced variate.
%6.       Fit a least squares line through your data (Christopher used a routine in Matlab)
%7.       The zero intercept (y = 0) is the modal value of your data
for col=1:size(hour_data,2)
    fit1=fit(sort_data(:,col),reduced_variate,'poly1');
    %plot(fit1,sort_data(:,col),reduced_variate)
    %refline(0,0)
    modal_maxima(1,col)=-fit1.p2/fit1.p1;
end


end

