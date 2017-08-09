function [ smoothed_data] = Running_mode( unsmoothed_data, window )
    %hopefully this will smoothe the data

    n=length(unsmoothed_data);
    n_wins=floor(n/window);
    for i=1:n_wins
        smoothed_data(((i-1)*window+1):i*window)=unsmoothed_data(((i-1)*window+1):i*window)-mode(unsmoothed_data(((i-1)*window+1):i*window));
    end
    remainder=n-(i*window+1);
    smoothed_data(i*window+1:i*window+1+remainder)=unsmoothed_data(i*window+1:i*window+1+remainder)-mode(unsmoothed_data(i*window+1:i*window+1+remainder));
    smoothed_data=smoothed_data';
end
