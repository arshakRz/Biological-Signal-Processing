%% BSP HW1 Arshak Rezvani 98106531
%% 0-1
clc; close all; clear all;
load(".\\Data_0\\EEG.mat");
fs_eeg = 512;
t = (1:length(eeg_signal)) * 1/fs_eeg;
figure
subplot(4, 1, 1);
plot(t, eeg_signal);
xlim([0 1]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The First Second of The Signal (Mean = ",num2str(mean(eeg_signal)),", STD = ", num2str(std(eeg_signal)),")");
title(tit, 'Interpreter','latex')
grid minor
eeg_signal_norm = (eeg_signal-mean(eeg_signal))/std(eeg_signal);
subplot(4, 1, 2);
plot(t, eeg_signal_norm);
xlim([0 1]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The First Second of The Normalized Signal (Mean = ",num2str(mean(eeg_signal_norm)),", STD = ", num2str(std(eeg_signal_norm)),")");
title(tit, 'Interpreter','latex')
grid minor
bpFilt = designfilt('bandpassfir','FilterOrder',50, ...
         'CutoffFrequency1',14,'CutoffFrequency2',30, ...
         'SampleRate',fs_eeg);
eeg_signal_bp_filter = filter(bpFilt, eeg_signal_norm);
subplot(4, 1, 3);
plot(t, eeg_signal_bp_filter);
xlim([0 1]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The First Second of The Bandpass Signal Using filter");
title(tit, 'Interpreter','latex')
grid minor
eeg_signal_bp_filtfilt = filtfilt(bpFilt, double(eeg_signal_norm));
subplot(4, 1, 4);
plot(t, eeg_signal_bp_filtfilt);
xlim([0 1]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The First Second of The Bandpass Signal Using filtfilt");
title(tit, 'Interpreter','latex')
grid minor
fvtool(bpFilt)
%% 0-2
clc; close all; clear all;
load(".\\Data_0\\ECG.mat");
fs_ecg = 128;
ecg_signal = ecg_singal;
clear ecg_singal;
t = (1:length(ecg_signal)) * 1/fs_ecg;
figure
subplot(5,1,1);
plot(t, ecg_signal);
t_max = 7;
t_min = 1.5;
xlim([t_min t_max]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The Signal between t = 1.5 to t = 7");
title(tit, 'Interpreter','latex')
grid minor
bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
         'HalfPowerFrequency1',5,'HalfPowerFrequency2',12, ...
         'SampleRate',fs_ecg);
ecg_signal_bp = filtfilt(bpFilt, ecg_signal);
subplot(5,1,2);
plot(t, ecg_signal_bp);
xlim([t_min t_max]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The First Ten Seconds of The Bandpass Signal");
title(tit, 'Interpreter','latex')
grid minor
z = tf('z');
H = (-z^(-2)-2*z^(-1)+2*z^1+z^2);
Num = H.Numerator;
ecg_signal_diff = filtfilt(Num{:}/10,1,ecg_signal_bp);
subplot(5,1,3);
plot(t, ecg_signal_diff);
xlim([t_min t_max]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The First Ten Seconds of The Diff Signal");
title(tit, 'Interpreter','latex')
grid minor
ecg_signal_diff_pow2 = ecg_signal_diff.^2;
subplot(5,1,4);
plot(t, ecg_signal_diff_pow2);
xlim([0 t_max]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The First Ten Seconds of The Diff Signal To The Power of Two");
title(tit, 'Interpreter','latex')
grid minor
ecg_signal_MA = movmean(ecg_signal_diff_pow2,20);
subplot(5,1,5);
plot(t, ecg_signal_MA);
xlim([t_min t_max]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("The First Ten Seconds of The MA Signal");
title(tit, 'Interpreter','latex')
grid minor
hold on 
Rs = zeros(size(ecg_signal_MA));
Rs(ecg_signal_MA>max(ecg_signal_MA)*0.8) = 100;
plot(t, Rs);
ylim([min(ecg_signal_MA) max(ecg_signal_MA)])
subplot(5, 1, 1);
hold on;
plot(t, Rs*max(ecg_signal));
ylim([min(ecg_signal) max(ecg_signal)])
%% 1-1
clc; close all; clear all; 
fs = 100;
t = 1/fs:1/fs:10;
x = sin(2*pi*(5*t+2*t.^2));
figure
subplot(2,1,1);
plot(t, x);
xlim([1 1.2]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("Signal Between t = 1 to t = 1.2");
title(tit, 'Interpreter','latex')
grid minor
subplot(2,1,2);
plot(t, x);
xlim([4 4.2]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
tit = strcat("Signal Between t = 4 to t = 4.2");
title(tit, 'Interpreter','latex')
grid minor
df = (fs)/(length(t)-1);
f = -fs/2:df:fs/2;
fx = fftshift(fft(x))*1/fs;
figure
subplot(2,1,1);
plot(f, abs(fx));
xlabel('f(Hz)', 'Interpreter','latex')
ylabel('Magnitude','Interpreter','latex')
tit = strcat("Fourier Transform Magnitude");
title(tit, 'Interpreter','latex')
grid minor
subplot(2,1,2);
plot(f, angle(fx));
xlabel('f(Hz)', 'Interpreter','latex')
ylabel('Phase','Interpreter','latex')
tit = strcat("Fourier Transform Phase");
title(tit, 'Interpreter','latex')
grid minor
noise = randn(size(x));
x_noisy = x + noise;
fx_noisy = fftshift(fft(x_noisy))*1/fs;
figure
subplot(2,1,1);
plot(f, abs(fx_noisy));
xlabel('f(Hz)', 'Interpreter','latex')
ylabel('Magnitude','Interpreter','latex')
tit = strcat("Fourier Transform Magnitude");
title(tit, 'Interpreter','latex')
grid minor
subplot(2,1,2);
plot(f, angle(fx_noisy));
xlabel('f(Hz)', 'Interpreter','latex')
ylabel('Phase','Interpreter','latex')
tit = strcat("Fourier Transform Phase");
title(tit, 'Interpreter','latex')
grid minor
energy = zeros(4, 39);
for i = 1:39
    ts = (i:(50+(i-1)*25));
    energy(1, i) = bandpower(x_noisy(ts),fs,[0.5 4]);
    energy(2, i) = bandpower(x_noisy(ts),fs,[4 8]);
    energy(3, i) = bandpower(x_noisy(ts),fs,[8 13]);
    energy(4, i) = bandpower(x_noisy(ts),fs,[13 30]);
end
figure
subplot(4, 1, 1)
plot(1:39, energy(1,:))
xlabel('Window Number', 'Interpreter','latex')
ylabel('Energy','Interpreter','latex')
tit = strcat("$\delta$");
title(tit, 'Interpreter','latex')
grid minor
xlim([1 39])
subplot(4, 1, 2)
plot(1:39, energy(2,:))
xlabel('Window Number', 'Interpreter','latex')
ylabel('Energy','Interpreter','latex')
tit = strcat("$\theta$");
title(tit, 'Interpreter','latex')
grid minor
xlim([1 39])
subplot(4, 1, 3)
plot(1:39, energy(3,:))
xlabel('Window Number', 'Interpreter','latex')
ylabel('Energy','Interpreter','latex')
tit = strcat("$\beta$");
title(tit, 'Interpreter','latex')
grid minor
xlim([1 39])
subplot(4, 1, 4)
plot(1:39, energy(4,:))
xlabel('Window Number', 'Interpreter','latex')
ylabel('Energy','Interpreter','latex')
tit = strcat("$\alpha$");
title(tit, 'Interpreter','latex')
grid minor
xlim([1 39])
tit = strcat("Energy in Different Windows for Bandpower:");
sgtitle(tit, 'Interpreter','latex')
fs = 512;
load(".\\Data_1\\EEG_rest.mat");
energy = zeros(4, 39);
for i = 1:39
    ts = (i:(50+(i-1)*25));
    energy(1, i) = bandpower(EEG_rest(ts),fs,[0.5 4]);
    energy(2, i) = bandpower(EEG_rest(ts),fs,[4 8]);
    energy(3, i) = bandpower(EEG_rest(ts),fs,[8 13]);
    energy(4, i) = bandpower(EEG_rest(ts),fs,[13 30]);
end
figure
subplot(4, 1, 1)
plot(1:39, energy(1,:))
xlabel('Window Number', 'Interpreter','latex')
ylabel('Energy','Interpreter','latex')
tit = strcat("$\delta$");
title(tit, 'Interpreter','latex')
grid minor
xlim([1 39])
subplot(4, 1, 2)
plot(1:39, energy(2,:))
xlabel('Window Number', 'Interpreter','latex')
ylabel('Energy','Interpreter','latex')
tit = strcat("$\theta$");
title(tit, 'Interpreter','latex')
grid minor
xlim([1 39])
subplot(4, 1, 3)
plot(1:39, energy(3,:))
xlabel('Window Number', 'Interpreter','latex')
ylabel('Energy','Interpreter','latex')
tit = strcat("$\beta$");
title(tit, 'Interpreter','latex')
grid minor
xlim([1 39])
subplot(4, 1, 4)
plot(1:39, energy(4,:))
xlabel('Window Number', 'Interpreter','latex')
ylabel('Energy','Interpreter','latex')
tit = strcat("$\alpha$");
title(tit, 'Interpreter','latex')
grid minor
xlim([1 39])
tit = strcat("Energy in Different Windows for Bandpower:");
sgtitle(tit, 'Interpreter','latex')
%% 1-2 
clc; close all; clear all;
load(".\\Data_0\\ECG.mat");
fs_ecg = 128;
ecg_signal = ecg_singal;
clear ecg_singal;
t = (1:length(ecg_signal)) * 1/fs_ecg;
Rs = R_peaks(ecg_signal, fs_ecg);
figure
plot(t, ecg_signal);
t_max = 60;
t_min = 50;
xlim([t_min t_max]);
xlabel('t(S)', 'Interpreter','latex')
ylabel('Signal','Interpreter','latex')
num_R = sum(diff(Rs)==1);
HR = num_R/t(end) * 60;
tit = strcat("The ECG Signal and R peaks (HR = 71.4 bpm)");
title(tit, 'Interpreter','latex')
grid minor
hold on;
plot(t, Rs*max(ecg_signal));
ylim([min(ecg_signal) max(ecg_signal)])
[pxx, f] = pwelch(ecg_signal, [], [], [], fs_ecg);
figure
plot(f, pxx)
tit = strcat("PSD");
title(tit, 'Interpreter','latex')
xlabel('Frequency (Hz)', 'Interpreter','latex')
grid minor
%% 2
clc; close all; clear all;
x = rand([1 7000]);
%x = randn([1 7000]);
%a = rand([1 7000]);
%b = rand([1 7000]);
%x = cos(2*pi*b) .* sqrt(-2*log(1-a));
fs = 512;
t = (1:length(x)) * 1/fs;
figure
tit = "Samples Generated Using rand";
sgtitle(tit, 'Interpreter','latex')
subplot(5, 1, 1)
plot(t, x)
xlim([t(1) t(end)])
xlabel('Time( S)', 'Interpreter','latex')
tit = strcat("Signal (Sampling Frequency = 512 Hz)");
title(tit, 'Interpreter','latex')
grid minor
subplot(5, 1, 2)
histogram(x)
tit = strcat("Histogram");
title(tit, 'Interpreter','latex')
xlabel('Range', 'Interpreter','latex')
[pxx, f] = pwelch(x, [], [], [], fs, 'centered');
grid minor
subplot(5, 1, 3)
plot(f, pxx)
tit = strcat("PSD");
title(tit, 'Interpreter','latex')
xlabel('Frequency (Hz)', 'Interpreter','latex')
xlim([-1 1])
grid minor
[r, lags] = xcorr(x,['biased']);
subplot(5, 1, 4)
plot(lags/fs, r)
xlim([lags(1)/fs lags(end)/fs])
xlabel('Lag (S)', 'Interpreter','latex')
tit = strcat("Autocorrelation (Sampling Frequency = 512 Hz)");
title(tit, 'Interpreter','latex')
grid minor
fr = fftshift(fft(r))*1/fs;
subplot(5, 1, 5)
df = (fs)/(length(lags)-1);
f = -fs/2:df:fs/2;
plot(f, abs(fr))
tit = strcat("Fourier Transform of Autocorrelation");
title(tit, 'Interpreter','latex')
xlabel('Frequency (Hz)', 'Interpreter','latex')
xlim([-1 1])
grid minor
%% 3-1
clc; close all; clear all; 
NN = [100 500 1000 5000];
figure
sgtitle('Autocorrelation with Different Number of Lags', 'Interpreter','latex')
for p = 1:length(NN)
    subplot(4,1,p)
    x = randn(1, 10000);
    z = tf('z');
    a = 0.5;
    H = 1/(1-a*z^(-1));
    Num = H.Numerator;
    Den = H.Denominator;
    y = filter(Num{:},Den{:},x);
    N = NN(p);
    [r_est, lags] = xcorr(y,N,'unbiased');
    plot(lags, r_est);
    fs = 1;
    nn = -(N-1):1/fs:(N-1);
    r = (1-a.^(2*(N-abs(nn))))/(1-a^2).*a.^abs(nn);
    hold on
    plot(nn, r);
    tit = strcat("N = ", num2str(N));
    title(tit, 'Interpreter','latex')
    xlabel('Sample', 'Interpreter','latex')
    grid minor
    legend('xcorr','Analytic')
    xlim([-100 100])
end
%% 3-2-1
clc; close all; clear all;
load(".\\Data_0\\EEG.mat");
fs = 512;
%t = (1:length(eeg_signal)) * 1/fs_eeg;
[r_est, lags] = xcorr(eeg_signal,'unbiased');
plot(lags/fs, r_est);
tit = strcat("Autocorrelation");
title(tit, 'Interpreter','latex')
xlabel('Lag (S)', 'Interpreter','latex')
grid minor
xlim([-2 2])
%% 3-2-2
clc; close all; clear all;
load(".\\Data_2\\ECGPCG.mat");
t = (1:length(ECG)) * 1/fs;
figure
subplot(2, 1, 1);
plot(t, ECG)
xlim([3.5 8])
grid on
xlabel('Time (S)', 'Interpreter','latex')
tit = strcat("ECG");
title(tit, 'Interpreter','latex')
subplot(2, 1, 2);
plot(t, PCG)
xlim([3.5 8])
grid on
xlabel('Time (S)', 'Interpreter','latex')
tit = strcat("PCG");
title(tit, 'Interpreter','latex')

figure
subplot(2,1,1)
[r_est, lags] = xcorr(ECG,'unbiased');
plot(lags, r_est);
tit = strcat("Autocorrelation of ECG");
title(tit, 'Interpreter','latex')
xlabel('Lag (S)', 'Interpreter','latex')
%xlim([-5000 5000])
grid minor
subplot(2,1,2)
[r_est, lags] = xcorr(PCG,'unbiased');
plot(lags, r_est);
tit = strcat("Autocorrelation of PCG");
title(tit, 'Interpreter','latex')
xlabel('Lag (S)', 'Interpreter','latex')
%xlim([-5000 5000])
grid minor
%% 4
clc; close all; clear all; 
load(".\\Data_3\\EEG_p300.mat");
x = EEG_p300.data(31,:);
% target 8
idx = find(EEG_p300.markers_seq==8 & EEG_p300.markers_target == 1);
p = 0; 
for i = 1:length(idx)
    p = p + x((idx(i)-103):(idx(i)+512));
end
p = p/i;
fs = 512;
t = linspace(-0.2,1,length(p));
figure
subplot(4,1,1)
plot(t, p);
xlim([-0.2 1])
xlabel('Time (S)', 'Interpreter','latex')
tit = strcat("Average Response for Target Character 8");
title(tit, 'Interpreter','latex')
grid minor
% target 3
idx = find(EEG_p300.markers_seq==3 & EEG_p300.markers_target == 1);
p = 0;
for i = 1:length(idx)
    p = p + x((idx(i)-103):(idx(i)+512));
end
p = p/i;
subplot(4,1,2)
plot(t, p);
xlim([-0.2 1])
xlabel('Time (S)', 'Interpreter','latex')
tit = strcat("Average Response for Target Character 3");
title(tit, 'Interpreter','latex')
grid minor
% non-target 3
idx = find(EEG_p300.markers_seq==3 & EEG_p300.markers_target ==2);
p = 0;
for i = 1:length(idx)
    p = p + x((idx(i)-103):(idx(i)+512));
end
p = p/i;
subplot(4,1,3)
plot(t, p);
xlim([-0.2 1])
xlabel('Time (S)', 'Interpreter','latex')
tit = strcat("Average Response for Non-target Character 3");
title(tit, 'Interpreter','latex')
grid minor
% shift 
idx = find(EEG_p300.markers_seq==8 & EEG_p300.markers_target == 1);
ref = x((idx(1)-103):(idx(1)+512));
p = ref;
for i = 2:length(idx)
    sidx = 0;
    R = 0;
    for j = -100:100
        sig = x((idx(i)-103-j):(idx(i)+512-j));
        Rn = corrcoef(ref,sig);
        Rn = Rn(1, 2);
        if Rn>R
            sidx = j;
            R = Rn;
        end
    end
    p = p + x((idx(i)-103-sidx):(idx(i)+512-sidx));
end
p = p/length(idx);
subplot(4,1,4)
plot(t, p);
xlim([-0.2 1])
xlabel('Time (S)', 'Interpreter','latex')
tit = strcat("Average Shift Optimized Response for Target Character 8");
title(tit, 'Interpreter','latex')
grid minor