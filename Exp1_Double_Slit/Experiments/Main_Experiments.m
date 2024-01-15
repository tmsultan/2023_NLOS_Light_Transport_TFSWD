close all; clear all

%% From Raw Data

%% Plotting Parameters
LW = 1.5;
XY_Text = 14;
Title_Text = 16;
Number_Text = 12;

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

%% Create Folder
fname = '../Figures/Experimental_Results/';
if ~exist(fname)
    mkdir(fname)
end

LegendEntry = {'With Diffuser', 'No Diffuser'}


%% Load and Trim Data
num1 = 37; % Occluded: Options: 18 - 28, except 22
num2 = 38; % UnOccluded


% Trimming because other reflections were detected
% start = 19150;
% stop = 19500;

start = 18900;
stop = 19500;

% Timing Axis
Time_Axis = 4*(1:stop - start + 1);


% start = 1;
% stop = 65536;

occdata = load(['data_Run_', num2str(num1),'.mat'], 'data3');
occ = occdata.data3(:,:,start:stop);
save([fname, 'With_Occluder_Trimmed.mat'], 'occ')
occ_int = squeeze(sum(occ, 1))./(max(occ(:)));


unoccdata = load(['data_Run_', num2str(num2),'.mat'], 'data3');
unocc = unoccdata.data3(:,:,start:stop);
save([fname, 'Without_Occluder_Trimmed.mat'], 'unocc')
unocc_int = squeeze(sum(unocc, 1))./(max(unocc(:)));

occ_int_sum = sum(occ_int(:));
unocc_int_sum = sum(unocc_int(:));

% Save for Noise Characteristics
save([fname,'Photon_Count'], "unocc_int_sum", "occ_int_sum" );


occ = squeeze(occ);
unocc = squeeze(unocc);

%% Visualize Dataset

figure; imagesc(Time_Axis, X, squeeze(occ));
ylabel('Detector Plane (cm)','FontSize', XY_Text)
xlabel('Time (ps)','FontSize', XY_Text)
hold on; title(LegendEntry{1}, 'FontSize', XY_Text);
colorbar; ax = gca;
ax.FontSize = Number_Text; 
hold off
saveas(gcf, [fname, 'Dataset_Occluded'], 'png');
savefig(gcf, [fname, 'Dataset_Occluded']);


figure; imagesc(Time_Axis, X, squeeze(unocc));
ylabel('Detector Plane (cm)','FontSize', XY_Text)
xlabel('Time (ps)','FontSize', XY_Text)
hold on; title(LegendEntry{2}, 'FontSize', XY_Text);
colorbar; ax = gca;
ax.FontSize = Number_Text; 
hold off
saveas(gcf, [fname, 'Dataset_UnOccluded'], 'png');
savefig(gcf, [fname, 'Dataset_UnOccluded']);

%% Zero Pad in Time Domain to reduce ringing
N = stop - start + 1;
% padsize = (1681 - N)/2;
% occ_int = padarray(occ_int,[0, padsize]);
% unocc_int = padarray(unocc_int,[0, padsize]);
% N = size(occ_int,2)

%% Isolate Pfield Component


tmp1 = 1/(4*10^-12)/N;
%PF = (6*10^9)/tmp1;
L5 = 4;
F5 = (3*10^8)/(L5*10^-2);
PF = (3*10^8)/(L5*10^-2)/tmp1;

L6 = 6;
F6 = (3*10^8)/(L6*10^-2);
PF2 = (3*10^8)/(L6*10^-2)/tmp1;

%PF = (7.5*10^9)/tmp1;


% Sx = size(occdata.data3,2);
% ticx= round([0:5:Sx]);



Y = fft((occ_int), [], 2);
Y2 = fft((unocc_int), [], 2);

% Y = Y./max(Y, [], 1);
% Y2 = Y2./max(Y2, [], 1);

Norm_Const = max(max(abs(Y), [], 1), max(abs(Y2), [], 1))

Y = Y./Norm_Const;
Y2 = Y2./Norm_Const;

figure
plot(X, abs(Y(:, 1)), 'LineWidth', LW);
hold on; 
plot(X, abs(Y2(:,1)), 'LineWidth', LW);
title('DC Component', 'FontSize', Title_Text)
legend(LegendEntry, 'location', 'best')
xlabel('Detector Plane (cm)', 'FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
% set(gca,'xtick', 1:4:Sx, 'xTickLabel',ticx)
ax = gca;
ax.FontSize = Number_Text; 
hold off
saveas(gcf, [fname, 'Zero_Freq'], 'png');
savefig(gcf, [fname, 'Zero_Freq']);



figure
plot(X, abs(Y(:, round(PF) - 1)),'LineWidth', LW)
hold on; 
plot(X, abs(Y2(:,round(PF)  - 1)), 'LineWidth', LW)
MSG = [num2str(F5/10^9), ' GHz/' num2str(L5) ' cm Wavelength'];
title(MSG, 'FontSize', Title_Text)
xlabel('Detector Plane (cm)','FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
legend(LegendEntry, 'location', 'best')
ax = gca;
ax.FontSize = Number_Text; 
hold off
saveas(gcf, [fname, '4_cm'], 'png');
savefig(gcf, [fname, '4_cm']);


 
figure
plot(X, abs(Y(:, round(PF2) - 1)),'LineWidth', LW)
hold on; 
plot(X, abs(Y2(:,round(PF2)  - 1)),'LineWidth', LW)
MSG2 = [num2str(F6/10^9), ' GHz/' num2str(L6) ' cm Wavelength'];
title(MSG2, 'FontSize', XY_Text);
xlabel('Detector Plane (cm)','FontSize', XY_Text)
ylabel('Normalized Intensity', 'FontSize', XY_Text)
legend(LegendEntry, 'location', 'best')
ax = gca;
ax.FontSize = Number_Text; 
hold off
saveas(gcf, [fname, '6_cm'], 'png');
savefig(gcf, [fname, '6_cm']);

N_Exp = N;
X_Exp = X;

save('../Results_Mat_Files/Experimental_Results_Occluded', 'occ', 'occ_int', 'Y', 'X_Exp', 'N_Exp')
save('../Results_Mat_Files/Experimental_Results_UnOccluded', 'unocc', 'unocc_int', 'Y2', 'X_Exp', 'N_Exp')
