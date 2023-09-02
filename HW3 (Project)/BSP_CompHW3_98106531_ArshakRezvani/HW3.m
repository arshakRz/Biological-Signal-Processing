%% BSP Computer HW3 98106531 Arshak Rezvani
%% Q1
clc; close all; clear all; 
addpath("Q1\");
load('ECG.mat');
ecg_signal = ecg_singal;
clear ecg_singal
fs = 128;
t = (1:length(ecg_signal))/fs;
% A
figure
subplot(3,1,1)
plot(t, ecg_signal)
title("The Original ECG","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 8])
sig_noisy = ecg_signal + 20*cos(2*pi*50*(t));
subplot(3,1,2)
plot(t, sig_noisy)
title("The Original ECG Contaminated With 50Hz Line Noise","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 8])
M = 5;
w = randn(1, M)';
d = sig_noisy;
x = 20*rand*cos(2*pi*50*(t)+rand*2*pi)';
s_hat = zeros(size(x));
n_hat = zeros(size(x));
for i = M:length(x)
    w = w + 0.00001*(d(i)-w'*flip(x((i-M+1):i)))*flip(x((i-M+1):i));
    s_hat(i) = w'*flip(x((i-M+1):i));
    n_hat(i) = d(i)-s_hat(i);
end
subplot(3,1,3)
plot(t, n_hat)
title("The Filtered ECG","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 8])
% B
ecg_signal_rep = repmat(ecg_signal(337:462), 1, 20);
t = (1:length(ecg_signal_rep))/fs;
figure
subplot(3,1,1)
plot(t, ecg_signal_rep)
title("The Original Repeated ECG","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 12])
sig_noisy = ecg_signal_rep + 20*randn(size(ecg_signal_rep));
subplot(3,1,2)
plot(t, sig_noisy)
title("The Original Repeated ECG Contaminated With White Gaussian Noise","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 12])
M = 5;
w = randn(1, M)';
d = sig_noisy;
x = zeros(size(d));
shift = 125;
x(1:(end-shift)) = d(shift+1:end);
x = x';
s_hat = zeros(size(x));
n_hat = zeros(size(x));
for i = M:length(x)
    w = w + 0.00001*(d(i)-w'*flip(x((i-M+1):i)))*flip(x((i-M+1):i));
    s_hat(i) = w'*flip(x((i-M+1):i));
    n_hat(i) = d(i)-s_hat(i);
end
subplot(3,1,3)
plot(t, s_hat)
title("The Filtered ECG","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 12])
% C
t = (1:length(ecg_signal))/fs;
figure
subplot(3,1,1)
plot(t, ecg_signal)
title("The Original ECG","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 12])
sig_noisy = ecg_signal + 20*randn(size(ecg_signal));
subplot(3,1,2)
plot(t, sig_noisy)
title("The Original ECG Contaminated With White Gaussian Noise","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 12])
M = 10;
w = randn(1, M)';
d = sig_noisy;
x = zeros(size(d));
shift = 116;
x(1:(end-shift)) = d(shift+1:end);
x = x';
s_hat = zeros(size(x));
n_hat = zeros(size(x));
for i = M:length(x)
    w = w + 0.00000001*(d(i)-w'*flip(x((i-M+1):i)))*flip(x((i-M+1):i));
    s_hat(i) = w'*flip(x((i-M+1):i));
    n_hat(i) = d(i)-s_hat(i);
end
subplot(3,1,3)
plot(t, s_hat)
title("The Filtered ECG","Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
xlim([0 12])
%% Q2
clc; close all; clear all; 
addpath("Q2\");
load("Normal.mat")
load("PVC.mat")
load("PVC1.mat")
load("tinv.mat")

figure
subplot(4,3,1)
plot(ecg_signal)
xlim([1 length(ecg_signal)])
title("The Normal Cycle","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")

[~, ecg_signal_mn] = rceps(ecg_signal);
subplot(4,3,2)
plot(ecg_signal_mn)
xlim([1 length(ecg_signal_mn)])
title("The Normal Cycle Min Phase Part","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")

[~, nd] = cceps(ecg_signal);
ecg_signal_mx = icceps(cceps(ecg_signal)-cceps(ecg_signal_mn), nd);
subplot(4,3,3)
plot(ecg_signal_mx)
xlim([1 length(ecg_signal_mx)])
title("The Normal Cycle Max Phase Part","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")

subplot(4,3,4)
plot(PVC)
xlim([1 length(PVC)])
title("The PVC Cycle","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
[~, PVC_mn] = rceps(PVC);
subplot(4,3,5)
plot(PVC_mn)
xlim([1 length(PVC_mn)])
title("The PVC Cycle Min Phase Part","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
[~, nd] = cceps(PVC);
PVC_mx = icceps(cceps(PVC)-cceps(PVC_mn), nd);
subplot(4,3,6)
plot(PVC_mx)
xlim([1 length(PVC_mx)])
title("The PVC Cycle Max Phase Part","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")

subplot(4,3,7)
plot(PVC1)
xlim([1 length(PVC1)])
title("The PVC1 Cycle","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
[~, PVC1_mn] = rceps(PVC1);
subplot(4,3,8)
plot(PVC1_mn)
xlim([1 length(PVC1_mn)])
title("The PVC1 Cycle Min Phase Part","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
[~, nd] = cceps(PVC1);
PVC1_mx = icceps(cceps(PVC1)-cceps(PVC1_mn), nd);
subplot(4,3,9)
plot(PVC1_mx)
xlim([1 length(PVC1_mx)])
title("The PVC1 Cycle Max Phase Part","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
subplot(4,3,10)
plot(tinv)
xlim([1 length(tinv)])
title("The tinv Cycle","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
[~, tinv_mn] = rceps(tinv);
subplot(4,3,11)
plot(tinv_mn)
xlim([1 length(tinv_mn)])
title("The tinv Cycle Min Phase Part","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
[~, nd] = cceps(tinv);
tinv_mx = icceps(cceps(tinv)-cceps(tinv_mn), nd);
subplot(4,3,12)
plot(tinv_mx)
xlim([1 length(tinv_mx)])
title("The tinv Cycle Max Phase Part","Interpreter","latex")
xlabel("Sample", "Interpreter","latex")
ylabel("Amp", "Interpreter","latex")
%% Q3
clc; close all; clear all; 
addpath("Q3\");
load('train_recording.mat');
load('train_annotations.mat');
load('test_recording.mat');
load('test_annotations.mat');
PCGCellArray = train_recordings;
annotationsArray = train_annotations;
numberOfStates = 4;
numPCGs = length(PCGCellArray);
Fs = 1000;
% A matrix of the values from each state in each of the PCG recordings:
state_observation_values = cell(numPCGs,numberOfStates);
num_points = 0;
for PCGi = 1:length(PCGCellArray)
    num_points = num_points + length(PCGCellArray{PCGi})/20;
end
features_tot = zeros(4, num_points)';
states_tot = zeros(1, num_points)';
ind = 1; 
for PCGi = 1:length(PCGCellArray)
    PCG_audio = PCGCellArray{PCGi};
    S1_locations = annotationsArray{PCGi,1};
    S2_locations = annotationsArray{PCGi,2};    
    [PCG_Features, featuresFs] = getSpringerPCGFeatures(PCG_audio, Fs);
    PCG_states = labelPCGStates(PCG_Features(:,1),S1_locations, S2_locations, featuresFs, 0);
    features_tot(ind:(ind+length(PCG_states)-1), :) = PCG_Features;
    states_tot(ind:(ind+length(PCG_states))-1) = PCG_states;
    ind = ind + length(PCG_states) - 1;
end

PCG_audio = test_recordings{1};
S1_locations = test_annotations{1};
S2_locations = test_annotations{2};
[features_tot_test, featuresFs] = getSpringerPCGFeatures(PCG_audio, Fs);
states_tot_test = labelPCGStates(features_tot_test(:,1),S1_locations, S2_locations, featuresFs, 0);

sig = interp(PCG_Features(:,1), 20)';
states = zeros(size(sig));
for i = 1:length(PCG_states)
    states((i-1)*20+1:i*20) = PCG_states(i);
end
t = (1:length(sig))/Fs;
z = zeros(size(t));
col = states;  % This is the color, vary with x in this case.
surface([t;t],[sig;sig],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',0.6);
xlim([0 5])
colormap('jet')
grid minor
xlabel('Time (s)','Interpreter','latex')
ylabel('Amp','Interpreter','latex')
title("Signal With Annotations",'Interpreter','latex')
hold on 
qw{1} = plot(nan);
qw{2} = plot(nan);
qw{3} = plot(nan);
qw{4} = plot(nan); % You can add an extra element too
legend([qw{:}], {'State 1','State 2','State 3', 'State 4'});

features_selector = zeros(4, 2);

test_labels_pred = zeros(size(features_tot_test));

figure
for k = 1:4
    labels = states_tot;
    labels(labels==k) = 0;
    labels(labels~=0) = 1;
    labels = double(not(labels));
    score = 0; 
    for i = 1:3
        for j = (i+1):4
            if score < fisher_score_nd(features_tot(:, [i j]), labels)
                score = fisher_score_nd(features_tot(:, [i j]), labels);
                features_selector(k, :) = [i j];
            end
        end
    end 
    Mdl = fitcnb(features_tot(:, features_selector(k, :)), labels,'Prior','uniform','DistributionNames','normal');
    [test_labels_pred(:, k), post, ~] = predict(Mdl, features_tot_test(:, features_selector(k, :)));
    post = test_labels_pred(:, k) .* post;
    test_labels_pred(:, k) = max(post,[],2);
    subplot(2,2,k)
    X = features_tot(:, features_selector(k, :));
    gscatter(X(:,1), X(:,2), labels);
    tit = strcat("The Best Pair Of Features For Class ", num2str(k), " Vs. Others (Features ", num2str(features_selector(k, 1)), " and ", num2str(features_selector(k, 2)),")");
    title(tit, 'Interpreter','latex')
    grid minor
    h = gca;
    cxlim = h.XLim;
    cylim = h.YLim;
    hold on
    Params = cell2mat(Mdl.DistributionParameters); 
    Mu = Params(2*(1:2)-1,1:2); % Extract the means
    Sigma = zeros(2,2,2);
    c = ['r-','c-'];
    for j = 1:2
        Sigma(:,:,j) = diag(Params(2*j,:)).^2; % Create diagonal covariance matrix
        xlim = Mu(j,1) + 4*[-1 1]*sqrt(Sigma(1,1,j));
        ylim = Mu(j,2) + 4*[-1 1]*sqrt(Sigma(2,2,j));
        f = @(x,y) arrayfun(@(x0,y0) mvnpdf([x0 y0],Mu(j,:),Sigma(:,:,j)),x,y);
        fcontour(f,[xlim ylim],c(j)) % Draw contours for the multivariate normal distributions 
    end
    h.XLim = [-2 4];
    h.YLim = [-2 4];
    legend('0','1','','','Location','northeast')
    hold off
end
clear xlim ylim
[~, test_labels_pred] = max(test_labels_pred,[],2);
pred = zeros(4, length(test_labels_pred));
pred(1, test_labels_pred==1) = 1;
pred(2, test_labels_pred==2) = 1;
pred(3, test_labels_pred==3) = 1;
pred(4, test_labels_pred==4) = 1;
actual = zeros(4, length(states_tot_test));
actual(1, states_tot_test==1) = 1;
actual(2, states_tot_test==2) = 1;
actual(3, states_tot_test==3) = 1;
actual(4, states_tot_test==4) = 1;
figure
plotconfusion(actual,pred,'Test')

sig = test_recordings{1}';
states = zeros(size(sig));
for i = 1:length(test_labels_pred)
    states((i-1)*20+1:i*20) = test_labels_pred(i);
end
t = (1:length(sig))/Fs;
z = zeros(size(t));
col = states;  % This is the color, vary with x in this case.
figure
surface([t;t],[sig;sig],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',0.6);
xlim([0 5])
colormap('jet')
grid minor
xlabel('Time (s)','Interpreter','latex')
ylabel('Amp','Interpreter','latex')
title("Test Signal With Predicted Annotations",'Interpreter','latex')
hold on 
%{
sig = interp(features_tot_test(:,1), 20)';
surface([t;t],[sig;sig],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',0.6);
%}
qw{1} = plot(nan);
qw{2} = plot(nan);
qw{3} = plot(nan);
qw{4} = plot(nan); % You can add an extra element too
legend([qw{:}], {'State 1','State 2','State 3', 'State 4'});
%% Q4
clc; close all; clear all; 
addpath("Q3\");
load('train_recording.mat');
load('train_annotations.mat');
load('test_recording.mat');
load('test_annotations.mat');
PCGCellArray = train_recordings;
annotationsArray = train_annotations;
numberOfStates = 4;
numPCGs = length(PCGCellArray);
Fs = 1000;
% A matrix of the values from each state in each of the PCG recordings:
state_observation_values = cell(numPCGs,numberOfStates);
num_points = 0;
for PCGi = 1:length(PCGCellArray)
    num_points = num_points + length(PCGCellArray{PCGi})/20;
end
features_tot = zeros(4, num_points)';
states_tot = zeros(1, num_points)';
ind = 1; 
for PCGi = 1:length(PCGCellArray)
    PCG_audio = PCGCellArray{PCGi};
    S1_locations = annotationsArray{PCGi,1};
    S2_locations = annotationsArray{PCGi,2};    
    [PCG_Features, featuresFs] = getSpringerPCGFeatures(PCG_audio, Fs);
    PCG_states = labelPCGStates(PCG_Features(:,1),S1_locations, S2_locations, featuresFs, 0);
    features_tot(ind:(ind+length(PCG_states)-1), :) = PCG_Features;
    states_tot(ind:(ind+length(PCG_states))-1) = PCG_states;
    ind = ind + length(PCG_states) - 1;
end
PCG_audio = test_recordings{1};
S1_locations = test_annotations{1};
S2_locations = test_annotations{2};
[features_tot_test, featuresFs] = getSpringerPCGFeatures(PCG_audio, Fs);
states_tot_test = labelPCGStates(features_tot_test(:,1),S1_locations, S2_locations, featuresFs, 0);
trans = zeros(numberOfStates);
features_tot = features_tot(1:10568,:);
states_tot = states_tot(1:10568);
diff_states_tot = diff(states_tot);
diff_states_tot(diff_states_tot~=0) = 1;
diff_states_tot = [0; diff_states_tot];
states_tot_4 = zeros(size(states_tot));
states_tot_4(states_tot==4) = 1;
states_tot_3 = zeros(size(states_tot));
states_tot_3(states_tot==3) = 1;
states_tot_2 = zeros(size(states_tot));
states_tot_2(states_tot==2) = 1;
states_tot_1 = zeros(size(states_tot));
states_tot_1(states_tot==1) = 1;
trans(4, 4) = (length(find(states_tot==4))-1*sum(and(diff_states_tot, states_tot_4)))/length(find(states_tot==4));
trans(3, 3) = (length(find(states_tot==3))-1*sum(and(diff_states_tot, states_tot_3)))/length(find(states_tot==3));
trans(2, 2) = (length(find(states_tot==2))-1*sum(and(diff_states_tot, states_tot_2)))/length(find(states_tot==2));
trans(1, 1) = (length(find(states_tot==1))-1*sum(and(diff_states_tot, states_tot_1)))/length(find(states_tot==1));
trans(4, 1) = 1-trans(4, 4);
trans(3, 2) = 1-trans(3, 3);
trans(2, 4) = 1-trans(2, 2);
trans(1, 3) = 1-trans(1, 1);
B = mnrfit(features_tot, states_tot);
emis = mnrval(B, features_tot)';
emis_test = mnrval(B, features_tot_test)';
estimatedStates = hmmviterbi(1:length(states_tot), trans, emis);
%estimatedStates = myViterbi(1:length(states_tot), trans, emis, [0.25 0.25 0.25 0.25]);
estimatedStates_test = hmmviterbi(1:length(states_tot_test), trans, emis_test);
pred = zeros(4, length(estimatedStates));
pred(1, estimatedStates==1) = 1;
pred(2, estimatedStates==2) = 1;
pred(3, estimatedStates==3) = 1;
pred(4, estimatedStates==4) = 1;
actual = zeros(4, length(states_tot));
actual(1, states_tot==1) = 1;
actual(2, states_tot==2) = 1;
actual(3, states_tot==3) = 1;
actual(4, states_tot==4) = 1;
figure
plotconfusion(actual,pred,'Train')

pred = zeros(4, length(estimatedStates_test));
pred(1, estimatedStates_test==1) = 1;
pred(2, estimatedStates_test==2) = 1;
pred(3, estimatedStates_test==3) = 1;
pred(4, estimatedStates_test==4) = 1;
actual = zeros(4, length(states_tot_test));
actual(1, states_tot_test==1) = 1;
actual(2, states_tot_test==2) = 1;
actual(3, states_tot_test==3) = 1;
actual(4, states_tot_test==4) = 1;
figure
plotconfusion(actual,pred,'Test')
sig = interp(features_tot_test(:,1), 20)';
states = zeros(size(sig));
for i = 1:length(states_tot_test)
    states((i-1)*20+1:i*20) = estimatedStates_test(i);
end
t = (1:length(sig))/Fs;
z = zeros(size(t));
col = states;  % This is the color, vary with x in this case.
figure
surface([t;t],[sig;sig],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',0.6);
xlim([0 10])
colormap('jet')
grid minor
xlabel('Time (s)','Interpreter','latex')
ylabel('Amp','Interpreter','latex')
title("Signal With Annotations Test",'Interpreter','latex')
hold on 
qw{1} = plot(nan);
qw{2} = plot(nan);
qw{3} = plot(nan);
qw{4} = plot(nan); % You can add an extra element too
legend([qw{:}], {'State 1','State 2','State 3', 'State 4'});

sig = interp(features_tot(:,1), 20)';
states = zeros(size(sig));
for i = 1:length(states_tot)
    states((i-1)*20+1:i*20) = estimatedStates(i);
end
t = (1:length(sig))/Fs;
z = zeros(size(t));
col = states;  % This is the color, vary with x in this case.
figure
surface([t;t],[sig;sig],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',0.6);
xlim([0 10])
colormap('jet')
grid minor
xlabel('Time (s)','Interpreter','latex')
ylabel('Amp','Interpreter','latex')
title("Signal With Annotations Train",'Interpreter','latex')
hold on 
qw{1} = plot(nan);
qw{2} = plot(nan);
qw{3} = plot(nan);
qw{4} = plot(nan); % You can add an extra element too
legend([qw{:}], {'State 1','State 2','State 3', 'State 4'});