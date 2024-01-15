clear all; clc; close all


%% Plotting Parameters
LW = 1.5;
XY_Text = 14;
Title_Text = 16;
Number_Text = 12;

% Horizontal
Nx = 66; % Number of Sampling of points in Horizontal Dimension
Lx = 0.01; %m Total Length of Sampling (60cm)
%Sx = (Lx)/(Nx-1); %m  Space between Samples


% % Vertical Plane
% % Horizontal
Ny = 19; % Number of Sampling of points in Horizontal Dimension
Ly = 58.2/100; %m Total Length of Sampling (60cm)
% Sy = (Ly)/(Ny-1); %m  Space between Samples

%% Experimental Parameters
SLT_Pln = 1.16;
SlitWidth = 2/100;
SlitSpacing = 50/100;

OccPln = SLT_Pln + 1.18;
OccWidth = 0.105/2;
Shift = [-0.012, 0, 0];

DetPln = OccPln + 0.77;

X = linspace(-65/2, 65/2, 66)+2.5;

%% Simulation Parameters
voxelSize = 0.005; % m

% 1 == Include Correction Factors; 0 == Don't Include
Correction_Factor = 1; 

% Add Poisson Noise
NOISE = 0; 

%% Create Folder

fname = ['../Figures/Results_TFSWD/'];
if ~exist(fname)
    mkdir(fname)
end

LegendEntry = {'With Diffuser', 'No Diffuser'}



% 1D Simulator

LP_Wall = [0, 0, 0];

% Ignore LP --> Shift it
%LP_Origin = [-1.5, 0, 1];
%D1 = dist(LP_Wall, LP_Wall)

% Coordinates_Target_Plane
% Slit 1
leftbottomfront = [-SlitSpacing/2 - SlitWidth,  -0.45, SLT_Pln] + Shift;
righttopback = [-SlitSpacing/2, 0.45, SLT_Pln] + Shift;
[Slit_1, volumeSize] = makegrid(leftbottomfront, righttopback, voxelSize);


leftbottomfront = [SlitSpacing/2, -0.45, SLT_Pln] + Shift;
righttopback = [SlitSpacing/2 + SlitWidth, 0.45, SLT_Pln] + Shift;
[Slit_2, volumeSize] = makegrid(leftbottomfront, righttopback, voxelSize);




%% Occluder
% leftbottomfront = [-0.055, -0.45, OccPln];
% righttopback = [0.055, 0.45, OccPln];
% [Occluder, volumeSize] = makegrid(leftbottomfront, righttopback, voxelSize);



%% Detection Plane
leftbottomfront = [X(1)/100, 0, DetPln];
righttopback = [X(end)/100, 0, DetPln];
[DP, volumeSize] = makegrid(leftbottomfront, righttopback, voxelSize);


X_Sim = [-32.5/100:voxelSize:32.5/100].*100;
% X_Sim = [leftbottomfront(1):voxelSize:righttopback(1)].*100;

%% Intialize Timebins
Slits = [Slit_1; Slit_2];
c = 3*10^8;
BinRes = 60*10^-12;
Range = round(4/(c*4*10^-12));
dataset.data = zeros(size(DP,1), 1, Range);
dataset.deltat = (c*4*10^-12);
dataset.t = Range;

c = physconst('LightSpeed'); % Loading speed of light
numpixel = BinRes/dataset.deltat*c; % Number of pixels at Full Width Half Max
sigma = numpixel/2*sqrt(2*log(2)); % Gaussian std%% Intialize Timebins

%% Core Algorithm
j = 1;
ct = 0;
% for loop with ray Intersection
Normal_Vector = [0,0,1];
for i = 1:size(Slits,1)
    Point1 = Slits(i, :);
    Dist1 = sqrt(sum((Point1 - LP_Wall).^2));
    cosfactor_1 = dot([Point1 - LP_Wall], Normal_Vector)./(norm([Point1 - LP_Wall]));
    for k = 1:size(DP,1)
        
        Point2 = DP(k, :);
        cosfactor_2 = dot([Point2-Point1], Normal_Vector)./(norm([Point2-Point1]));
        
        %%
        %         % Check if it intersects occluder
        %         m = (Point2(1) - Point1(1))/(Point2(3) - Point1(3));
        %         c = Point1(1) - m.*Point1(3);
        %
        %         % Check Intersection
        %         x3 = m.*OccPln+c;
        %         %if ((abs(x3) < 0.05) || (x3 >= 0.05))
        %
        %         if (abs(x3) < 0)
        %             continue
        %         end
        %%
        
        Dist2 = sqrt(sum((Point1 - Point2).^2));
        D = Dist1 + Dist2;
        
        % Create the target
        index = floor(double(D)/double(dataset.deltat));
        if (index > 0) && (index < dataset.t)
            
            % Correction Factors
            if Correction_Factor == 1
                albedo = 1; distfactor = (1./(Dist1*Dist2))^2; 
                cosfactor = cosfactor_1*cosfactor_2;
            else
                albedo = 1; distfactor = 1; cosfactor = 1;
            end
            
            dataset.data(k, j, :) = squeeze(dataset.data(k, j, :)) + ...
                squeeze(albedo*distfactor*cosfactor*exp(-1/2.*(((1:dataset.t)-index)/sigma).^2))';
        end
        
        ct = ct + 1;
        
        
    end
end

% Visualize Dataset
figure;
plot(squeeze(dataset.data(5,:,:)))
figure;
imagesc(squeeze(dataset.data))
hold on; title('Without Occluder'); hold off

dataset1 = dataset;



%% Rerun Occluded Case
dataset_Occ = dataset;
dataset_Occ.data = zeros(size(DP,1), 1, Range);

% Define Occluder Plane
n = [0, 0, 1]; % Plane Normal
P0 = [0 , 0, OccPln]; % Point in Occ Plane

j = 1;
ct = 0;
% for loop with ray Intersection
for i = 1:size(Slits,1)
    Point1 = Slits(i, :);
    Dist1 = sqrt(sum((Point1 - LP_Wall).^2));
    cosfactor_1 = dot([Point1 - LP_Wall], Normal_Vector)./(norm([Point1 - LP_Wall]));

    for k = 1:size(DP,1)
        
        Point2 = DP(k, :);
        cosfactor_2 = dot([Point2-Point1], Normal_Vector)./(norm([Point2-Point1]));
        
        
        %% Calculate whether ray is blocked
        % Wikipedia: https://en.wikipedia.org/wiki ...
        % /Line%E2%80%93plane_intersection
        
        % Calculate Direction Vector
        L = Point2 - Point1;
        
        
        % Calculate Scaling Factor
        if (dot(L, n) == 0)
            display('Error: Occluder Plane is Parallel to Optical Axis')
            break
        end
        
        d_SF = dot((P0 - Point1), n)./(dot(L, n));
        
        % Calculate intersection point
        IntPt_OccPln = Point1 + L.*d_SF;
        
        %         % Check if it intersects occluder
        %         m = (Point2(1) - Point1(1))/(Point2(3) - Point1(3));
        %         c = Point1(1) - m.*Point1(3);
        
        %
        %         % Check Intersection
        %         x3 = m.*OccPln+c;
        %         %if ((abs(x3) < 0.05) || (x3 >= 0.05))
        %
        if ((abs(IntPt_OccPln(1)) < OccWidth) && (abs(IntPt_OccPln(2)) < 0.45))
            continue
        end
        
        Dist2 = sqrt(sum((Point1 - Point2).^2));
        D = Dist1 + Dist2;
        
        % Create the target
        index = floor(double(D)/double(dataset_Occ.deltat));
        if (index > 0) && (index < dataset_Occ.t)
            
            % Correction Factors
            if Correction_Factor == 1
                albedo = 1; distfactor = (1./(Dist1*Dist2))^2; 
                cosfactor = cosfactor_1*cosfactor_2;
            else
                albedo = 1; distfactor = 1; cosfactor = 1;
            end
            
            dataset_Occ.data(k, j, :) = squeeze(dataset_Occ.data(k, j, :)) + ...
                squeeze(albedo*distfactor*cosfactor*exp(-1/2.*(((1:dataset_Occ.t)-index)/sigma).^2))';
        end
        
        ct = ct + 1;
        
        
    end
end

% Visualize Dataset
figure;
plot(squeeze(dataset_Occ.data(5,:,:)))
figure;
imagesc(squeeze(dataset_Occ.data))
hold on; title('With Occluder'); hold off

%% Add Noise

%% Adding Poisson Noise
if NOISE == 1
    start = 1;
    stop = 3000;
    load('Photon_Count.mat')
    dcr = 0.001;
    uniformNoise = 0;
    ambientNoise = 0;
    dataset_Occ.data = noisyDataset(dataset_Occ.data, occ_int_sum, dcr, uniformNoise, ambientNoise)
    dataset1.data = noisyDataset(dataset1.data, unocc_int_sum, dcr, uniformNoise, ambientNoise)
end

%datasetNC.data = ones(1, 1000)/10;

% Max_Dataset = max(dataset.data(:));
% 
% % Scale
% dataset_noisy= dataset.data./Max_Dataset;
% J = imnoise(dataset_noisy,'poisson');
% 
% Max_Dataset2 = 217;
% 
% % Noisy Dataset
% J1 = round(J.*Max_Dataset2);
% dataset1.data = J1;
% %dataset.data = J1;
% 


%% Pfield Convolution
start = 3000-601+1;
stop = 3000;

start = 1650;
stop = 3330;




% start = 1;
% stop = Range;

% Timing Axis
Time_Axis = 4*(1:stop - start + 1);

occ = dataset_Occ.data(:,:,start:stop);
occ_int = squeeze(sum(occ, 2))./(max(occ(:)));
save([fname, 'With_Occluder_Trimmed_Simulated.mat'], 'occ')
occ = squeeze(occ);


unocc = dataset1.data(:,:,start:stop);
unocc_int = squeeze(sum(unocc, 2))./(max(unocc(:)));
save([fname, 'Without_Occluder_Trimmed_Simulated.mat'], 'unocc')
unocc = squeeze(unocc);




%% Visualization
figure; imagesc(Time_Axis, X_Sim, squeeze(occ));
ylabel('Detector Plane (cm)','FontSize', XY_Text)
xlabel('Time (ps)','FontSize', XY_Text)
hold on; title(LegendEntry{1}, 'FontSize', XY_Text);
colorbar; ax = gca;
ax.FontSize = Number_Text; 
hold off
saveas(gcf, [fname, 'Dataset_Occluded'], 'png');


figure; imagesc(Time_Axis, X_Sim, squeeze(unocc));
ylabel('Detector Plane (cm)','FontSize', XY_Text)
xlabel('Time (ps)','FontSize', XY_Text)
hold on; title(LegendEntry{2}, 'FontSize', XY_Text);
colorbar; ax = gca;
ax.FontSize = Number_Text; 
hold off
saveas(gcf, [fname, 'Dataset_UnOccluded'], 'png');





%% Isolate Pfield Component

N = stop - start + 1;
tmp1 = 1/(4*10^-12)/N;
%PF = (6*10^9)/tmp1;
L5 = 4;
F5 = (3*10^8)/(L5*10^-2);
PF = (3*10^8)/(L5*10^-2)/tmp1;

L6 = 6;
F6 = (3*10^8)/(L6*10^-2);
PF2 = (3*10^8)/(L6*10^-2)/tmp1;

%PF = (7.5*10^9)/tmp1;


%Sx = size(occdata.data3,2);
%ticx= round([0:5:Sx]);



Y = fft((occ_int), [], 2); Y = Y./max(Y, [], 1);
Y2 = fft((unocc_int), [], 2); Y2 = Y2./max(Y2, [], 1);

Y_Sim = Y; Y2_Sim = Y2;



figure; imagesc(Time_Axis, X_Sim, squeeze(abs(Y2)));
ylabel('Detector Plane (cm)','FontSize', XY_Text)
xlabel('Time (ps)','FontSize', XY_Text)
hold on; title(LegendEntry{1}, 'FontSize', XY_Text);
colorbar; ax = gca;
ax.FontSize = Number_Text; 
hold off
saveas(gcf, [fname, 'Dataset_Occluded'], 'png');


figure
plot(X_Sim, abs(Y(:, 1)), 'LineWidth', LW);
hold on; 
plot(X_Sim, abs(Y2(:,1)), 'LineWidth', LW);
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
plot(X_Sim, abs(Y(1:end, round(PF) - 1)),'LineWidth', LW)
hold on; 
plot(X_Sim, abs(Y2(1:end,round(PF)  - 1)), 'LineWidth', LW)
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
plot(X_Sim, abs(Y(1:end, round(PF2) - 1)),'LineWidth', LW)
hold on; 
plot(X_Sim, abs(Y2(1:end,round(PF2)  - 1)),'LineWidth', LW)
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



save('../Results_Mat_Files/Simulation_Results_TFSWD_Occluded', 'occ', 'occ_int', 'Y_Sim', 'X_Sim', 'N')
save('../Results_Mat_Files/Simulation_Results_TFSWD_UnOccluded', 'unocc', 'unocc_int', 'Y2_Sim', 'X_Sim', 'N')

%% Unoccluded case is same as Phasor Field
save('../Results_Mat_Files/Simulation_Results_PF_UnOccluded', 'unocc', 'unocc_int', 'Y2_Sim', 'X_Sim', 'N')