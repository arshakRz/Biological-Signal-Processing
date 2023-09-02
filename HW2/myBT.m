function [psd_x] = myBT(x, w, M)
    if w == "rectwin"
        x_window = x;
    elseif w == "kaiser"
        x_window = x.*kaiser(length(x), 38);
    else 
        x_window = x.*w;
    end
    R = xcorr(x_window)/M;
    psd_x = fftshift(fft(R));
end