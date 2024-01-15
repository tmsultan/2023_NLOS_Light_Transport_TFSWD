R = 29*2.54/2
OD = 88;
% F = 24*2.54;
% ID = 1/((1/F) - (1/ID))

ID = 172.2

F = OD*ID/(OD+ID)
% ID1 = 80.5 % (To Diffuser)
% ID2 = ID - ID1 % (Diffuser to Detector)
% ID2 = 80.7
% ID1 = ID- ID2

% F = 26*2.54


% R = 29*2.54/2
R = R*cos(atan(30/90))

N = stop - start + 1;
tmp1 = 1/(4*10^-12)/N;
L5 = 4;
F5 = (3*10^8)/(L5*10^-2);
PF = (3*10^8)/(L5*10^-2)/tmp1;
%PF = (6*10^9)/tmp1;
%PF = (7.5*10^9)/tmp1;



% Spot Size
% D_Spot = 2*1.22*(L5/100)*F/(R*2)

% Spot Size
D_Spot = 1.22*(L5/100)*ID/(R)


% Spot Size
D_Spot = 1.22*(L5/100)*ID/(R*2)


% Optical Spot
D_Opt = (0.5/100)*ID/OD;

D_Spot_T = sqrt(D_Spot^2+D_Opt^2)

% % Width of Gaussian
D_Spot = 2.355*0.84*(L5/100)*ID/(R*2) * 1/2




% Y_fit_opt = abs( Y(1:end, 1))./max(abs( Y(1:end,  1)));
Y_fit = abs( Y(1:end, round(PF) - 1))./max(abs( Y(1:end, round(PF) - 1)));


% F_fit = (Y(1:end, round(PF) - 1)./max(abs( Y(1:end, round(PF) - 1))));
% y = lowpass(F_fit,0.1);
% figure; plot(abs(y))
% Y_fit = y./max(abs(y));

% fit_list = []
% for i = 1:N
% Edge = data(i, 3e4:4e4);
figure; plot(abs(Y_fit))
f = fit(X', abs(Y_fit),'gauss1');

% fit_list = [fit_list; f.a1, f.b1, f.c1];
% figure; plot(f); hold on; plot(Edge)
% end


%% SIGMA == D_Spot*100/2.355
% FWHM = 2.355*Sigma
y = gaussmf(X, [D_Spot*100/2.355,  f.b1]);

% 

MSG = [num2str(round(F5/10^9, 1)), ' GHz/' num2str(L5) ' cm Wavelength'];
figure;
plot(X, abs(Y_fit), 'r--', 'LineWidth', 2 ) %Phasor Field is defined to be the fourier transform of the intensity
hold on; 
plot(X,y,'k', 'LineWidth', 2)
hold off
title(MSG)
xlabel('Detector Plane (cm)')
ylabel('Normalized Intensity')
legend('Experimental', 'Theoretical', 'location', 'best')


% 
% y = lowpass(Y_fit_opt,0.01);
% figure; plot(abs(y))
% Y_fit_opt = y./max(abs(y));




MSG = [num2str(round(F5/10^9, 1)), ' GHz/' num2str(L5) ' cm Wavelength'];
figure;
plot(X, abs(Y_fit), 'r--', 'LineWidth', 2 ) %Phasor Field is defined to be the fourier transform of the intensity
hold on; 
% plot(X,y,'k', 'LineWidth', 2)
plot(s_x, I_Airy, 'k', 'LineWidth', 2)
hold off
title(MSG)
xlabel('Detector Plane (cm)')
ylabel('Normalized Intensity')
legend('Experimental', 'Theoretical', 'location', 'best')



MSG = ['Phasor Field (', num2str(L5),' cm)']
figure;
plot(X,Y_fit_opt, 'LineWidth', 2)
hold on; 
plot(X, abs(Y_fit), 'r--', 'LineWidth', 2)
hold off
title('')
xlabel('Detector Plane (cm)')
ylabel('Normalized Intensity')
title('Focused Pfield Spot')
legend('Optical Carrier (DC)', MSG, 'location', 'best')


