clear;
clc;
% for real human PPI networks
dataset = 'protein_HuRI';
mlset = 'HuRI_ml';
k = [1000];
[G, n, pc] = semi_nmf_bmc(dataset, mlset, k);
