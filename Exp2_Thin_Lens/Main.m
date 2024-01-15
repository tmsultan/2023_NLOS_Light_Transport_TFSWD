close all; clear all; clc
%% Load File
load("DataProcessed\data_NoDiffuser.mat", 'Data_ND')
load("DataProcessed\data_WithDiffuser.mat", 'Data_WD')

N = size(Data_WD, 2);
%% Input Parameters
lambda = 4e-2; % Phasor Field Wavelength (cm)

c = physconst('lightspeed');
ts = 4e-12; % Temporal Sampling
% N = 101; % Number of Samples
X = linspace(-63.8/2, 63.8/2, 101);
xs = X(2) - X(1); % Spatial Sampling
%% Plotting Parameters
LW = 1.5;
XY_Text = 14;
Title_Text = 16;
Number_Text = 12+4;



%% Create Folder
fname = 'Results_Paper/';
if ~exist(fname)
    mkdir(fname)
end

%% Isolate DC and PF = 4cm Component from Measurement
% No Diffuser
[Y_PF_Norm_ND_DC, ~, ~, ~, ~] = PF_Component(Data_ND, N, 0, ts); % DC
[Y_PF_Norm_ND, ~, ~, Freq, ~] = PF_Component(Data_ND, N, lambda, ts);

% With Diffuser
[Y_PF_Norm_WD_DC, ~, ~, ~, ~] = PF_Component(Data_WD, N, 0, ts); % DC
[Y_PF_Norm_WD, ~, ~, Freq, ~] = PF_Component(Data_WD, N, lambda, ts);


%% Compute Theoretical 
% Airy_Bessel_Compute
Lambda_Opt = 5.15e-7; xs = Lambda_Opt/4;
[I_Airy_Opt, s_x_Opt] = Airy_Bessel_Compute(xs, Lambda_Opt);
I_Airy_Opt = circshift(I_Airy_Opt, -round(0.638/xs)); % shift maximum to experimental max


Lambda = 0.04; xs = Lambda/4;
[I_Airy, s_x] = Airy_Bessel_Compute(xs, Lambda);
I_Airy = circshift(I_Airy, -round(0.638/xs)); % shift maximum to experimental max

close all
%% 4x4 Plot
% Figure 11
figure
tiledlayout(2,2);

nexttile

plot(X, abs(Y_PF_Norm_WD), 'r--x', 'LineWidth', LW)
hold on; 
% title('Theoretical vs Experimental')
xlabel('Detector Plane (cm)', 'FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
ax.FontSize = Number_Text; 
plot(s_x, I_Airy, 'k', 'LineWidth', LW);
set(gca,'Box','on');
hold off
% legend('Experimental', 'Theoretical', 'location', 'best')
% legend('Experimental', 'Theoretical',  'location', 'east')


nexttile
hold on; 

plot(X, abs(Y_PF_Norm_ND), 'r--x', 'LineWidth', LW)
hold on; 
% title('Theoretical vs Experimental')
xlabel('Detector Plane (cm)', 'FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
ax.FontSize = Number_Text; 
plot(s_x, I_Airy, 'k', 'LineWidth', LW);
set(gca,'Box','on');
hold off
% legend('Experimental', 'Theoretical', 'location', 'best')
% legend('Experimental', 'Theoretical',  'location', 'east')

nexttile
hold on; 

plot(X, abs(Y_PF_Norm_WD), 'r--x', 'LineWidth', LW)
hold on; 
% title('Theoretical vs Experimental')
xlabel('Detector Plane (cm)', 'FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
ax.FontSize = Number_Text; 
plot(s_x, I_Airy, 'k', 'LineWidth', LW);
set(gca,'Box','on');
hold off
% legend('Experimental', 'Theoretical', 'location', 'best')
% legend('Experimental', 'Theoretical',  'location', 'east')



nexttile
hold on; 
plot(s_x_Opt, I_Airy_Opt, 'k', 'LineWidth', LW);
% plot(s_x, I_Airy, 'k', 'LineWidth', LW);
plot(X, abs(Y_PF_Norm_ND), 'r--x', 'LineWidth', LW)
xlabel('Detector Plane (cm)','FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
% legend('Theoretical', 'Experimental', 'location', 'southeast')
ax = gca;
ax.FontSize = Number_Text; 
set(gca,'Box','on');
% plot(s_x, I_Airy, 'k', 'LineWidth', LW);
hold off
legend('Theoretical', 'Experimental', 'location', 'southeast')

% 

fig = gcf;
fig.Position(3) = 1.75*(fig.Position(3))
fig.Position(4) = 2*(fig.Position(4))

% add legend
Lgnd = legend('show');
% Lgnd.Position(1) = 0.03;% 0.454761906997079;
% Lgnd.Position(2) = 0.471031746031746;
Lgnd.FontSize = Number_Text*1.25;
Lgnd.Layout.Tile = 'South' 
% Lgnd.Layout.Tile = 'East' 


movegui(gcf, 'onscreen')

saveas(gcf,  [fname, '11_Combined_TL_4_cm'], 'png');
% savefig(gcf, [fname, '11_Combined_TL_4_cm'], 'compact');


%% Plot Zero Frequency
% Figure 12(a)
figure
plot(X, abs(Y_PF_Norm_WD_DC), 'LineWidth', LW);
hold on; 
% title('Focused P-Field Spot', 'FontSize', Title_Text)
xlabel('Detector Plane (cm)', 'FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
ax.FontSize = Number_Text; 
plot(X, abs(Y_PF_Norm_ND_DC) , 'LineWidth', LW)
hold off
% legend('Optical Carrier (DC)', 'Phasor Field (4 cm)', 'location', 'best')
legend('With Diffuser', 'No Diffuser', 'location', 'best')

fig = gcf;
fig.Position(3) = 1.75/2*(fig.Position(3))
fig.Position(4) = 1*(fig.Position(4))

saveas(gcf, [fname, '12_Focused_PField_Spot_DC'], 'png' )
savefig(gcf, [fname, '12_Focused_PField_Spot_DC'] )



%% Create 1x2 Plot: Multiple Wavelengths
% Figure 13
[Y_PF_Norm_ND_3, ~, ~, Freq3, ~] = PF_Component(Data_ND, N, 0.03, ts);

figure; 
subplot(1, 2, 2)
plot(X, abs(Y_PF_Norm_ND_DC), 'LineWidth', LW);
hold on; 
% sgtitle('Focused P-Field Spot', 'FontSize', Title_Text)
title('No Diffuser')

xlabel('Detector Plane (cm)', 'FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
ax.FontSize = Number_Text; 
plot(X, abs(Y_PF_Norm_ND_3), '--', 'LineWidth', LW)
plot(X, abs(Y_PF_Norm_ND), 'k--', 'LineWidth', LW)

hold off
% legend('DC', '3cm',  '4cm', 'location', 'best')
legend('DC', '3cm',  '4cm', 'location', 'east')
% legend('DC', '3cm',  '4cm', 'location', 'eastoutside')

% legend('Optical Carrier (DC)', 'Phasor Field (4 cm)',  'location', 'east')

[Y_PF_Norm_WD_5, ~, ~, Freq2, ~] = PF_Component(Data_WD, N, 0.05, ts);
[Y_PF_Norm_WD_3, ~, ~, Freq3, ~] = PF_Component(Data_WD, N, 0.03, ts);



subplot(1, 2, 1)

plot(X, abs(Y_PF_Norm_WD_DC), 'LineWidth', LW)
hold on; 
title('With Diffuser')
xlabel('Detector Plane (cm)', 'FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
ax.FontSize = Number_Text; 
% plot(X, abs(Y_PF_Norm_WD_5), 'LineWidth', LW)
plot(X, abs(Y_PF_Norm_WD_3), 'LineWidth', LW)
plot(X, abs(Y_PF_Norm_WD), 'k', 'LineWidth', LW)
hold off
% legend('Experimental', 'Theoretical', 'location', 'best')
legend('DC', '3cm',  '4cm', 'location', 'east')
% legend('DC', '3cm',  '4cm', 'location', 'eastoutside')

fig = gcf;
fig.Position(3) = 2*fig.Position(3)
fig.Position(4) = 1*fig.Position(4)
movegui(gcf, 'onscreen')


fname = 'Results_Paper/';
saveas(gcf, [fname, '13_PField_MW'], 'png' )
savefig(gcf, [fname, '13_PField_MW'], 'compact')
