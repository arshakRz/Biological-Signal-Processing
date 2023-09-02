function [Rs] = R_peaks(ecg_signal, fs)
    bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
         'HalfPowerFrequency1',5,'HalfPowerFrequency2',12, ...
         'SampleRate',fs);
    ecg_signal_bp = filtfilt(bpFilt, ecg_signal);
    z = tf('z');
    H = (-z^(-2)-2*z^(-1)+2*z^1+z^2);
    Num = H.Numerator;
    ecg_signal_diff = filtfilt(Num{:}/10,1,ecg_signal_bp);
    ecg_signal_diff_pow2 = ecg_signal_diff.^2;
    ecg_signal_MA = movmean(ecg_signal_diff_pow2,20);
    Rs = zeros(size(ecg_signal_MA));
    Rs(ecg_signal_MA>max(ecg_signal_MA)*0.5) = 1;
end