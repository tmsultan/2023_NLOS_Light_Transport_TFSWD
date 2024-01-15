%% Plotting Parameters
LW = 1.5;
XY_Text = 14;
Title_Text = 16;
Number_Text = 12;



%%
% Parameters
z = 2;            % [=] m:  Depth Location of Point Source
a = 1;            % [=] m:   Depth Location of Occluders
z_offset = 1;     % [=] m:   Depth offset of Point source to Occluders
D = 2;            % [=] m:   Relay wall Diameter
gamma = 50e-12;   % [=] s:   Temporal precison of ToF detector
c = 3e8;          % [=] m/s: Speed of Light
lambda0 = 5e-7;   % [=] m:   Wavelength of Optical 
ds = 1e-3;        % [=] m:   Spatial extent of light source
x = 0:0.001:5;    % [=] m:   Transverse Coordinate

%% Plot Blur Kerne;
% FWHM
f1 = 3^(1/6)*z;                                 % Intensity Fall-Off
f2 = 2*c*gamma*(sqrt(0.5*D^2+z^2)/D);           % Time of Flight
f3 = sqrt(2*lambda0*(z_offset/a*(a+z_offset))+(z_offset*ds/a)^2);    % Occlusions/Shadows

% Standard Deviation
f1 = f1/(2*sqrt(2*log(2)))
f2 = f2/(2*sqrt(2*log(2)))
f3 = f3/(2*sqrt(2*log(2)))

% Gaussain
T1 = exp(-x.^2./(2*f1^2));
T2 = exp(-x.^2./(2*f2^2));
T3 = exp(-x.^2./(2*f3^2));

 
% Log Plot - Blur Kernel
x = log10(x)
figure; plot(x, T1)
hold on; 
plot(x, T2)
plot(x, T3)
hold off
xlabel('log10(\rho)')
legend('Falloff', 'Timing', 'Shadows')


%% Fourier Domain (FD)
f = -60:0.01:60; % Frequency Range

% Standard Dev in Fourier Domain
ff1 = 1/(2*pi*f1);
ff2 = 1/(2*pi*f2);
ff3 = 1/(2*pi*f3);

% Gaussians in Fourier Domain
GG1 = exp(-2 * (pi*f*f1).^2);
GG2 = exp(-2 * (pi*f*f2).^2);
GG3 = exp(-2 * (pi*f*f3).^2).*(1./f);

% Plot
figure; plot(f, GG1,  'LineWidth', LW)
hold on; 
plot(f, GG2,  'LineWidth', LW)
plot(f, abs(GG3),  'LineWidth', LW)
ax = gca;
set(gca,'Box','on');
ax.FontSize = Number_Text; 
hold off
xlabel('f')
ylim([0,1])
legend('Fall off', 'Timing', 'Shadows')
% legend('\Delta x_n', '\Delta x_o', 'I(\rho)')


xlabel('k', 'FontSize', XY_Text)
ylabel('|F\{I(\rho)\}|', 'FontSize', XY_Text)
saveas(gcf, ['FourierAnalysis'], 'png')
