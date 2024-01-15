clc; clear all; close all


%% Load dataset without Diffuser 2

% filename = ['Results/2022_04_15_Results/NoDiffuser_101cmPlane_L17/Run_1/'] % No Diffuser, but peak is not sharp

load('Data\data_Raw_NoDiffuser.mat')
% Data_ND = 
%% Load dataset with Diffuser 2

% filename = ['Results/2022_04_15_Results/NoDiffuser_101cmPlane_L17/Run_1/'] % No Diffuser, but peak is not sharp


% filename = ['Results/2022_03_27_Results/Small_Mirror/Diffuser_L23/'] % Multiple peaks saptail noise
% filename = ['Results/2022_03_27_Results/Small_Mirror/Diffuser_L23_2/'] % Noise looks good, but dual peaks and mall area
% filename = ['Results/2022_03_27_Results/Small_Mirror/Diffuser_L23_5_Long/']
% filename = ['Results/2022_03_27_Results/Small_Mirror/Diffuser_L23_5/']
% filename = ['Results/2022_03_27_Results/Small_Mirror/Run_2/Diffuser_L23_6/']
% filename = ['Results/2022_04_15_Results/NoDiffuser_40cmPlane_L19/Run_2/'] % No Diffuser, but peak is not sharp
% filename = ['Results/2022_04_15_Results/NoDiffuser_40cmPlane_L17/Run_2/'] % No Diffuser, but peak is not sharp
% filename = ['Results/2022_04_15_Results/NoDiffuser_101cmPlane_L17/Run_1/'] % No Diffuser, but peak is not sharp
filename = ['Results/2022_04_15_Results/NoDiffuser_101cmPlane_L17_Moved1/Run_1/'] % No Diffuser, but peak is not sharp
% filename = ['Results/2022_04_15_Results/Diffuser_101cm_L25/Run_1/'] % No Diffuser, but peak is not sharp
% filename = ['Results/2022_04_15_Results/Diffuser_101cm_L25/Run_2/'] % No Diffuser, but peak is not sharp
% 
 filename = ['Results/2022_04_15_Results/Diffuser_101cm_L27/Run_1/'] %
% Best one we have 

% filename = ['Results/2022_04_15_Results/NoDiffuser_101cmPlane_L17/Run_2/'] % No Diffuser, but peak is not sharp
% filename = ['Results/2022_04_15_Results/NoDiffuser_101cmPlane_L17/Run_1/'] % No Diffuser, but peak is not sharp



N = 101
data = zeros(N, 65536);


for i = 1:N
    j = i-1;
    
    if j < 10
        MSG = [filename,'data   ' , num2str(j), '.dat'];
        data(i,:) = str2num(fileread(MSG));
    elseif j < 100
        MSG = [filename,'data  ' , num2str(j), '.dat'];
        data(i,:) = str2num(fileread(MSG));
    else
        MSG = [filename,'data ' , num2str(j), '.dat'];
        data(i,:) = str2num(fileread(MSG));
    end

    %data(i,) = str2num(fileread('data 200.dat'))

end


data_Raw = data;
imagesc(data_Raw)
save([filename, '\data_Raw.mat'], 'data_Raw')

% Process Data (find t0)

data_sum = sum(data, 1);
figure; 
plot(data_sum)


Manual_t0 = 26000;%31255% 30646;  %31255

% Manual_t0 = 30646;


start = Manual_t0;
stop = Manual_t0 + 1000*26;



data_crop = data(:, start:stop);
[M, I] = max(data_crop, [], 2);

X1 = 25.6%
X = linspace(0, 63.8, N);

%  X = 1:N
figure; plot(X, M)
title('Maximum Intensity')
xlabel('X (cm)')

% fit_list = []
% for i = 1:N
% Edge = data(i, 2.55e4:3e4);
% f = fit([1:length(Edge)]', Edge.','gauss1');
% fit_list = [fit_list; f.a1, f.b1, f.c1];
% % figure; plot(f); hold on; plot(Edge)
% end
% 
% figure; plot(X, fit_list(:, 3)); title('Gaussian Widths')
% xlabel('Lateral Distance (cm)')
% 
% figure; plot(X, fit_list(:, 2)); title('Gaussian Peak Location')
% xlabel('Lateral Distance (cm)')



% figure; 

% return
% % data2 = reshape(data, [19 19 length(data(1, :))]);
% % [M, I] = max(data, [], 2);
% 
% 
% % Check
% data3 = data;
% %data3(:, 1:2:end,:) = flipud(data2(:,1:2:end,:));
% 
% % Produce Fixed Data
% data3(:, 1:2:end,:) = flipud(data2(:,1:2:end,:));
% 
% M2 = reshape(M, [19 19])
% I2 = reshape(I, [19, 19])
% 
% % I2Copy = I2;
% %
% % figure; imagesc(I2); title('I2')
% %
% % I2(2:2:end,:) = flipud(I2(2:2:end,:))
% % I3 = I2Copy;
% 
% 
% I3 = I2;
% %I3(:, 2:2:end) = flipud(I3(:,2:2:end))
% I3(:, 1:2:end) = flipud(I3(:,1:2:end))
% figure; imagesc(I3); title('I3')





% %% Make Video
% timevec = round(linspace(11640, 11710,70))
% tic = 0:10:90;
% ticx= round([0:2:19].*5*70/90);
% ticy= round([0:2:19].*5*65/90);
% 
% 
% figure;
% for i = 1:length(timevec)
%     j = timevec(i);
%     %pause(0.5);
%     imagesc(data3(:,:,j))
%     set(gca,'xtick', 1:2:19, 'xTickLabel',ticx)
%     set(gca,'ytick', 1:2:19, 'yTickLabel',flipud(ticy')')
%     xlabel('cm'); ylabel('cm');
%     title(['time index: ' num2str(j)])
%     colorbar
% end
% 
% 
% figure(5); gif('myfileRun4_Fixed.gif', 'frame', figure(5));
% 
% hold on;
% set(gca,'xtick', 1:2:19, 'xTickLabel',ticx)
% set(gca,'ytick', 1:2:19, 'yTickLabel',(ticy))
% xlabel('cm'); ylabel('cm');
% for i = 1:length(timevec)
%     j = timevec(i);
%     %pause(0.5);
%     imagesc(data3(:,:,j))
%     title(['time index: ' num2str(j)])
%     colorbar
%     gif
% end
% 
% hold off;


%% Phasor Field
% num = 22;
% 
% start = t0;
% stop = 34428;
% 
% % start = 1;
% % stop = 65536;
% 
% occdata = load(['data_Run_', num2str(num),'.mat'], 'data3');
% occ = occdata.data3(:,:,start:stop);

%save('With_Occluder_Trimmed.mat', 'occ')


% occ_int = squeeze(sum(occ, 1));


N = stop - start + 1;
tmp1 = 1/(4*10^-12)/N;
L5 = 4;
F5 = (3*10^8)/(L5*10^-2);
PF = (3*10^8)/(L5*10^-2)/tmp1;
%PF = (6*10^9)/tmp1;
%PF = (7.5*10^9)/tmp1;



Y = fft((data_crop), [], 2);

figure
plot(X, abs(Y(:, 1)));
hold on; 
title('Optical - Zero Frequency')
xlabel(' Detector Plane (cm)')
hold off

%saveas(gcf, 'Zero_Freq', 'png');

% figure
% plot(abs(Y_int))
% hold on; 
% plot(abs(Y2(:,1)))
% hold off


figure
plot(X, abs(Y(1:end, round(PF) - 1)))
hold on;
MSG = [num2str(F5/10^9), ' GHz/' num2str(L5) ' cm Wavelength'];
%title('6 GHz/5 cm Wavelength')
title(MSG)
xlabel('Transverse Detector Plane (cm)')
hold off

%saveas(gcf, '5_cm', 'png');

return
Y2 =  abs(Y(1:end, round(PF) - 1))
Y2 = reshape(Y2, [7, 21])
Y2( 2:2:end,:) = flip(Y2( 2:2:end,:) )
figure; imagesc(Y2)