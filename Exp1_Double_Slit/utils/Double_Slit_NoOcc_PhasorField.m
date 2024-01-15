clear;
close all;


lambda=4e-2; %wavelength
c = physconst('LightSpeed'); 
k=2*pi/lambda; %wave number


%% Plotting Parameters
LW = 1.5;
XY_Text = 14;
Title_Text = 16;
Number_Text = 12;

% % Horizontal
% Nx = 66; % Number of Sampling of points in Horizontal Dimension
% Lx = 0.01; %m Total Length of Sampling (60cm)
% %Sx = (Lx)/(Nx-1); %m  Space between Samples
% 
% X = linspace(-60/2, 60/2, 66);
% 
% % % Vertical Plane
% % % Horizontal
% Ny = 19; % Number of Sampling of points in Horizontal Dimension
% Ly = 58.2/100; %m Total Length of Sampling (60cm)
% % Sy = (Ly)/(Ny-1); %m  Space between Samples

%% Create Folder
fname = 'Sent_July_22_2022_Ver2/';
if ~exist(fname)
    mkdir(fname)
end


%% Initialize with the E-Field at plane A
samples=10000; %samples in each plane
width=4; %width of the planes in the x dimension
dx=(2*width)/(samples-1); %dx for the numerical integrals

% xA goes from -2 mm to +2mm in increments of dx ~ 0.2 microns
xA=linspace(-width,width,samples);%x coordinate for the sample points in plane a
xA1 = xA - 0.055;
EA=zeros(1,samples)+1i*zeros(1,samples);
Spot_Size = 0.01 % 0.5 cm
EA(samples/2-round(Spot_Size/2/dx):samples/2+round(Spot_Size/2/dx))=1;%E-Field at plane A


xB=linspace(-width,width,samples);%x coordinate for the sample points in plane B
zB=1.16; %z coordinate for the sample points in plane B
EB=zeros(1,samples)+zeros(1,samples)*1i; %E-Field at plane B

%% Plot the E-Field
figure; 
subplot(1,2,1)
hold on;
plot(xA.*100,real(EA));
plot(xA.*100,imag(EA));
hold off;
title('E-Field at plane A'); 
xlabel('x position (cm)'); ylabel('Normalized E-Field'); legend('Real Part', 'Imaginary Part');
subplot(1,2,2)
plot(xA.*100,abs(EA.^2)); title('Intensity at plane A'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane A



%% Evaluate the propagation integral from plane A to plane B
disp('Propagating A -> B ...')

for l=1:samples
    r=sqrt((xB-xA(l)).^2+(zB-0)^2);%Distance between source and sample points
    EB=EB+EA(l).*(exp(1i*k*r))./r*dx; %integral (note the difference between 1 and l)
end

EB=EB./max(EB);%re-normalizing. Let's not worry about the coefficients

%% Plot the E-Field
figure; 
subplot(1,2,1)
hold on;
plot(xB.*100,real(EB));
plot(xB.*100,imag(EB));
hold off;
title('E-Field at plane B before mask'); 
xlabel('x position (cm)'); ylabel('Normalized E-Field'); legend('Real Part', 'Imaginary Part');
subplot(1,2,2)
plot(xB.*100,abs(EB.^2)); title('Intensity at plane B before mask'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B



%% Apply the mask at plane B
mask=zeros(1,samples);

Slit_Width = 2/100;
Slit_Spacing = 50/100;

% Slit 1
ind1 = find(min(abs(xB - Slit_Spacing/2)) == abs(xB-Slit_Spacing/2));
ind2 = find(min(abs(xB - Slit_Spacing/2 - Slit_Width)) == ...
    abs(xB - Slit_Spacing/2 - Slit_Width));

% Slit 2
ind3 = find(min(abs(xB + Slit_Spacing/2)) == abs(xB+Slit_Spacing/2));
ind4 = find(min(abs(xB + Slit_Spacing/2 + Slit_Width)) == ...
    abs(xB + Slit_Spacing/2 + Slit_Width));

mask(ind1:ind2)=1; % mask goes from -1.5 mm to +1.5 mm
mask(ind4:ind3)=1; % mask goes from -1.5 mm to +1.5 mm
figure; plot(xB,mask); title('Mask tansmission'); xlabel('x position (cm)'); 

EB_mask=EB.*mask; % Propagate for entire EB and then apply mask



%% Plot the E-Field
figure; 
subplot(1,2,1)
hold on;
plot(xB.*100,real(EB_mask));
plot(xB.*100,imag(EB_mask));
hold off;
title('E-Field at plane B after mask'); 
xlabel('x position (cm)'); ylabel('Normalized E-Field'); legend('Real Part', 'Imaginary Part');
subplot(1,2,2)
plot(xB.*100,abs(EB_mask.^2)); title('Intensity at plane B after mask'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B


%% Propagate from plane B to plane C
xC=linspace(-width,width,samples);%x coordinate for the sample points in plane C
% xC = linspace(-0.35, 0.35, 71);
samples1 = length(xC);

EC=zeros(1,samples1)+zeros(1,samples1)*1i;

% C is 10 cm away from B
zC=zB+1.18+0.77;%z locaton of destination plane (source plane is at z=0)
zC=zB+1.18+0.77;%z locaton of destination plane (source plane is at z=0)


disp('Propagating B -> C ...')

for l=1:samples
    r=sqrt((xC-xB(l)).^2+(zC-zB)^2);%Distance between source and sample points
    EC=EC+EB_mask(l).*exp(1i*k*r)./((zB-zC))*dx; %integral (note the difference between 1 and l)
end
% EC = EC./(1i.*lambda);
EC=EC./max(EC);%re-normalizing. Let's not worry about the coefficients

%% plot the result
figure; 
subplot(1,2,1)
hold on;
plot(xC.*100,real(EC));
plot(xC.*100,imag(EC));
hold off;
title('E-Field at plane C'); 
xlabel('x position (cm)'); ylabel('Normalized E-Field'); legend('Real Part', 'Imaginary Part');
subplot(1,2,2)
plot(xC.*100,abs(EC.^2)); title('Intensity at plane C'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B
% plot(xC.*100,abs(EC)); title('Intensity at plane C'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B
% plot(xC.*100,(EC.^2)); title('Intensity at plane C'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B
% plot(xC.*100,real(EC)); title('Intensity at plane C'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B

%% Plot all 
% Plot Field at A
figure; 
subplot(4,2,1)
hold on;
plot(xA.*100,real(EA));
plot(xA.*100,imag(EA));
hold off;
title('E-Field at plane A'); 
xlabel('x position (cm)'); ylabel('Normalized E-Field'); legend('Real Part', 'Imaginary Part');
subplot(4,2,2)
plot(xA.*100,abs(EA.^2)); title('Intensity at plane A'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane A


% Plot the E-Field at B before Mask
subplot(4,2,3)
hold on;
plot(xB.*100,real(EB));
plot(xB.*100,imag(EB));
hold off;
title('E-Field at plane B before mask'); 
xlabel('x position (cm)'); ylabel('Normalized E-Field'); legend('Real Part', 'Imaginary Part');
subplot(4,2,4)
plot(xB.*100,abs(EB.^2)); title('Intensity at plane B before mask'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B



% Plot the E-Field at B after mask

subplot(4,2,5)
hold on;
plot(xB.*100,real(EB_mask));
plot(xB.*100,imag(EB_mask));
hold off;
title('E-Field at plane B after mask'); 
xlabel('x position (cm)'); ylabel('Normalized E-Field'); legend('Real Part', 'Imaginary Part');
subplot(4,2,6)
plot(xB.*100,abs(EB_mask.^2)); title('Intensity at plane B after mask'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B



% plot the result
subplot(4,2,7)
hold on;
plot(xC.*100,real(EC));
plot(xC.*100,imag(EC));
hold off;
title('E-Field at plane C'); 
xlabel('x position (cm)'); ylabel('Normalized E-Field'); legend('Real Part', 'Imaginary Part');
subplot(4,2,8)
plot(xC.*100,abs(EC.^2)); title('Intensity at plane C'); xlabel('x position (cm)'); ylabel('Normalized Intensity') %plot intensity at plane B


save('../Results_Mat_Files/Simulation_Results_PF_UnOccluded', 'EC', 'xA', 'xA1', 'samples');


