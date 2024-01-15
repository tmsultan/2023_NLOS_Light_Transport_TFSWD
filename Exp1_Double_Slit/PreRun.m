clear all; clc; close all

%% Run this first to generate the mat files for the final figure

addpath("utils\")

%% Run Occluded Phasor Field RSD Model 
run('utils\Double_Slit_Occ_PhasorField_RSD.m')

%% Run Simulation (TFSWD)
run('utils\Double_Slit_TFSWD.m')


%% Generate Experimental Results
run(['Experiments\Main_Experiments.m'])