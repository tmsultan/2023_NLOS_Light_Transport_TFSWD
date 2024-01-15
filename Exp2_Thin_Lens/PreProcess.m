clc; clear all; close all


%% Trim Parameters (Crop relevent transient signal)
Manual_t0 = 26000;
start = Manual_t0;
stop = Manual_t0 + 1000*26;

ts = 4e-12; % Temporal Sampling
N = 101; % Number of Samples
X = linspace(-63.8/2, 63.8/2, 101);
xs = X(2) - X(1); % Spatial Sampling
%% Load dataset without Diffuser 2
% filename in original location
%   filename = ['Results/2022_04_15_Results/NoDiffuser_101cmPlane_L17/Run_1/'] 

load('DataRaw\data_Raw_NoDiffuser.mat')
Data_ND = data_Raw(:, start:stop); % ND [=] No Diffuser

save('DataProcessed\data_NoDiffuser.mat', 'Data_ND', 'data_Raw', 'ts', 'xs')

clear data_Raw Data_ND
%% Load dataset with Diffuser 2
% filename in original location
%   filename = ['Results/2022_04_15_Results/Diffuser_101cm_L27/Run_1/'] %

load('DataRaw\data_Raw_WithDiffuser.mat')
Data_WD = data_Raw(:, start:stop); % WD [=] With Diffuser

save('DataProcessed\data_WithDiffuser.mat', 'Data_WD', 'data_Raw', 'ts', 'xs')