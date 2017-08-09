function [ smoothing_line] = Running_mode1( unsmoothed_data, window )
    %hopefully this will smoothe the data

    n=length(unsmoothed_data);
    n_wins=floor(n/window);
    for i=1:n_wins
        smoothing(((i-1)*window+1):i*window)=mode(unsmoothed_data(((i-1)*window+1):i*window));
        %smoothed_data(((i-1)*window+1):i*window)=unsmoothed_data(((i-1)*window+1):i*window)-);
    end
    remainder=n-(i*window+1);
    smoothing(i*window+1:i*window+1+remainder)=mode(unsmoothed_data(i*window+1:i*window+1+remainder));
    smoothing_line=smoothing;
end
