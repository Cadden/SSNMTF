clear;
clc;

alpha = 1;
ratios = [0.01:0.01:0.1];
nRpt = 1;

type = 1; % 0 for non-overlapping networks; 1 for overlapping networks

% non-overlapping LFR networks (type = 0)

% dataset = 'network10002';%(u=0.6)
% dataset = 'network10003'; %(u=0.7)
% semi_nmf_LFR(dataset, type, alpha, ratios, nRpt);

% overlapping LFR networks (type = 1)

% dataset = 'network10006';%(u=0.6,on=200, om = 2)
dataset = 'network10007'; %(u=0.7, on=200,om = 2)
semi_nmf_LFR(dataset, type, alpha, ratios, nRpt);
