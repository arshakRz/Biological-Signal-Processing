function [psd_x] = myWelch(x, w, M, L)
    if w == "rectwin"
        x_window = x;
    elseif w == "kaiser"
        x_window = x.*kaiser(length(x), 38);
    else 
        x_window = x.*w;
    end
    X=reshape(x_window, 1, []);
    y=[];
    idx=L;
    i=0;
    while idx<=M
       yL= X(i+1: idx);
       y=[y; yL];
       i=i+L;
       idx=idx+L;
    end
    psd_x = abs(fftshift(fft(y, [], 2))).^2/L;
    psd_x = mean(psd_x, 1);

end