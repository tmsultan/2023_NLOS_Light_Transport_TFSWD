function [I_Airy, s_x] = Airy_Bessel_Compute(xs, Lambda)

%% Input
% xs: Sampling Rate  (m)
% Lambda: wavelength (m)

s_x = -32:xs:32;
ID = 172.2;         % Image Distance


%% Account for non-zero angle from optical axis
theta = atan(s_x./ID); % Radians
% lambda = 0.04;
k = 2*pi/Lambda;
R = 29*2.54/2;
R = R*cos(atan(30/90))/100;


I0 = 1;             % Normalized Intensity

x = k*R*sin(theta);
I_Airy = I0*(2*besselj(1,x)./(x)).^2;

% figure; plot(x, I_Airy)
figure; plot(s_x, I_Airy)

end






