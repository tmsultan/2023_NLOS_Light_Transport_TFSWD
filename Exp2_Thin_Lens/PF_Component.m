function [Y_PF_Norm, Y_PF, PF_Index, Freq, Lambda] = PF_Component(Data, N, Lambda, ts)

%%
% Input:
% Data: Raw Data
% N: Number of Temporal Samples
% c: Speed of Light (m/s)
% Lambda = Wavelength (m)


%%
c = physconst('lightspeed');
Freq = c/Lambda; 
% Omega = 2*pi*Freq;
fs = 1./(N*ts); % Sampling Frequency
PF_Index = Freq/fs; % Index corresponding to the PF frequency


Y = fft((Data), [], 2);

if Lambda == 0
    Y_PF = Y(:, 1);
    PF_Index = 1;
    Lambda = 0;
else
    Y_PF = Y(:, round(PF_Index)-1);
end

Y_PF_Norm = Y_PF./max(Y_PF);

end