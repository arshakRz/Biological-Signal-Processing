%% Q1
clc; close all; clear all; 
addpath("Q_1\");
x1 = randn(1, 1000);

f1 = 0.05;
f2 = 0.40;
f3 = 0.45;
n = 1:1000;
phi1 = rand*2*pi;
phi2 = rand*2*pi;
phi3 = rand*2*pi;
x2 = 2*cos(2*pi*f1*n+phi1) + 2*cos(2*pi*f2*n+phi2) + 2*cos(2*pi*f3*n+phi3);

A = [1 -1.5 1.4]; 
x3 = filter(1, A, randn(1,1000));

pxx = myPeriodogram(x2, "rectwin");
pxx1 = periodogram(x2)';
figure
subplot(2,1,1)
plot(pxx)
subplot(2,1,2)
plot(pxx1)