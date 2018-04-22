function [ smoothed_data modes_out] = Running_mode( unsmoothed_data, window )
    %hopefully this will smooth the data

    n=length(unsmoothed_data);
    n_wins=floor(n/window);
    for i=1:n_wins
        this_mode=mode(unsmoothed_data(((i-1)*window+1):i*window));
        smoothed_data(((i-1)*window+1):i*window)=unsmoothed_data(((i-1)*window+1):i*window)-this_mode;
        modes_out(((i-1)*window+1):i*window)=this_mode*ones(window,1);
    end
    remainder=n-(i*window+1);
    this_mode=mode(unsmoothed_data(i*window+1:i*window+1+remainder));
    smoothed_data(i*window+1:i*window+1+remainder)=unsmoothed_data(i*window+1:i*window+1+remainder)-this_mode;
    modes_out(i*window+1:i*window+1+remainder)=this_mode*ones(length(i*window+1:i*window+1+remainder),1);
    smoothed_data=smoothed_data';
end
