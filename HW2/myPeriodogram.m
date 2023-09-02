function [psd_x] = myPeriodogram(x, w) 
    if w == "rectwin"
        x_window = x;
    elseif w == "kaiser"
        x_window = x.*kaiser(length(x), 38);
    else 
        x_window = x.*w;
    end
    psd_x = abs(fftshift(fft(x_window))).^2/length(x); 
end