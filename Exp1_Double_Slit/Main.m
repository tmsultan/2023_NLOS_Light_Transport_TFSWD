clear all; close all

%% Plotting Parameters
LW = 1.5;
XY_Text = 14;
Title_Text = 16;
Number_Text = 16;

% Horizontal
Nx = 66; % Number of Sampling of points in Horizontal Dimension
Lx = 0.01; %m Total Length of Sampling (60cm)
%Sx = (Lx)/(Nx-1); %m  Space between Samples

X = linspace(-65/2, 65/2, 66);

% % Vertical Plane
% % Horizontal
Ny = 19; % Number of Sampling of points in Horizontal Dimension
Ly = 58.2/100; %m Total Length of Sampling (60cm)
% Sy = (Ly)/(Ny-1); %m  Space between Samples

% Phasor field Parameters
L5 = 4;
F5 = (3*10^8)/(L5*10^-2);


%% Create Folder
fname = 'Figures/Final/';
if ~exist(fname)
    mkdir(fname)
end


%% Load Experimental Results
load("Results_Mat_Files\Experimental_Results_Occluded.mat")
load("Results_Mat_Files\Experimental_Results_UnOccluded.mat")

tmp1 = 1/(4*10^-12)/N_Exp;
PF_Exp = (3*10^8)/(L5*10^-2)/tmp1;


%% Load Simulation using Phasor Field
% Occluded
load('Results_Mat_Files\Simulation_Results_PF_Occluded.mat', 'ED', 'xA', 'xA1', 'samples');

Shift = 0.055; %Shift axis to match the maximum value
mask = logical((xA >= -0.325 + Shift).*(xA <= 0.325 + Shift));
mask2 = logical((xA1 >= -0.325).*(xA1 <= 0.325));
xA1 = circshift(xA, round(5.5/100/(xA(2)-xA(1))));
% mask3 = logical((x_mask >= -0.325).*(x_mask <= 0.325));


% Load Simulation for UnOccluded Case
load('Results_Mat_Files\Simulation_Results_PF_UnOccluded.mat')
tmp1 = 1/(4*10^-12)/N;
PF_Sim = (3*10^8)/(L5*10^-2)/tmp1;

%% Load Simulation using TFSWD
load('Results_Mat_Files\Simulation_Results_TFSWD_Occluded.mat')
load('Results_Mat_Files\Simulation_Results_TFSWD_UnOccluded.mat')
tmp1 = 1/(4*10^-12)/N;
PF_Sim = (3*10^8)/(L5*10^-2)/tmp1;


%% 4x4 Plot
% Figure 6(b)
figure(1);
tiledlayout(2,2);


nexttile
% Phasor Feild 
hold on; 
plot(X_Sim, abs(Y2_Sim(1:end, round(PF_Sim) - 1)),'LineWidth', LW)
plot(X, abs(Y2(:,round(PF_Exp)  - 1)), '--x', 'LineWidth', LW)
xlabel('Detector Plane (cm)','FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
ax.FontSize = Number_Text; 
set(gca,'Box','on');
hold off

nexttile
hold on; 


plot((xA1(mask)).*100, abs(ED(mask).^2)./(abs(ED(samples/2 + 1).^2)),'LineWidth', LW)
plot(X, abs(Y(:,round(PF_Exp)  - 1)), '--x', 'LineWidth', LW)
xlabel('Detector Plane (cm)','FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
ax.FontSize = Number_Text; 
set(gca,'Box','on');
hold off

%
nexttile
hold on; 

plot(X_Sim, abs(Y2_Sim(1:end, round(PF_Sim) - 1)),'LineWidth', LW)
plot(X, abs(Y2(:,round(PF_Exp)  - 1)), '--x', 'LineWidth', LW)
xlabel('Detector Plane (cm)','FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
ax = gca;
set(gca,'Box','on');
ax.FontSize = Number_Text; 
hold off

nexttile
hold on; 
plot(X_Sim, abs(Y_Sim(:,round(PF_Sim)  - 1)),  'LineWidth', LW)
plot(X, abs(Y(:,round(PF_Exp)  - 1)), '--x', 'LineWidth', LW)
xlabel('Detector Plane (cm)','FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
legend('Theoretical', 'Experimental', 'location', 'southeast')

ax = gca;
ax.FontSize = Number_Text; 
set(gca,'Box','on');
hold off



%% Save Figure
figure(1)
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

saveas(gcf, [fname, 'Combined_4_cm'], 'png');
savefig(gcf, [fname, 'Combined_4_cm']);

%% Figure  6(a)
PF_Sim = 2; PF_Exp = 2;

figure
hold on; 
plot(X_Sim, abs(Y_Sim(:,round(PF_Sim)  - 1)),  'LineWidth', LW)
plot(X, abs(Y(:,round(PF_Exp)  - 1)), '--x', 'LineWidth', LW)
% MSG = ['Radio Wave ' , newline ,  num2str(F5/10^9), ' GHz/', num2str(L5*100), ' cm Wavelength'];
% MSG = [num2str(F5/10^9), ' GHz/', num2str(L5*100), ' cm Wavelength'];
MSG = ['With Occluder'];
MSG = ['Light Field'];
% title(MSG, 'FontSize', Title_Text)
% xlim([-32.5, 32.5])
xlabel('Detector Plane (cm)','FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
%legend('Theoretical', 'Experimental', 'location', 'southeast')
ax = gca;
ax.FontSize = Number_Text; 
set(gca,'Box','on');
hold off


fig = gcf;
fig.Position(3) = 1.75/2*(fig.Position(3))
% fig.Position(4) = 2*(fig.Position(4))

movegui(gcf, 'onscreen')
saveas(gcf, [fname, 'TFSWD_', num2str(PF_Sim-2), '_cm'], 'png');
savefig(gcf, [fname, 'TFSWD_', num2str(PF_Sim-2), '_cm']);